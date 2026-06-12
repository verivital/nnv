function property = load_vnnlib2(propertyFile)
    % property = load_vnnlib2(propertyFile)
    %
    % Parse a VNN-LIB **2.0** specification file (header `(vnnlib-version <2.0>)`)
    % for the SINGLE-NETWORK case and return the SAME contract as load_vnnlib (1.0):
    %
    %     property.lb    -> input lower-bound vector (single) OR cell{M} of vectors
    %     property.ub    -> input upper-bound vector (single) OR cell{M} of vectors
    %     property.prop  -> cell{N} of structs, each with field .Hg = HalfSpace array
    %                       (HalfSpace.G*y <= HalfSpace.g encodes the UNSAFE output
    %                        region; an array of length>1 is an OR of disjuncts,
    %                        exactly as load_vnnlib produces).
    %
    % SOUNDNESS (the VNN-COMP -150 rule): any 2.0 construct this parser cannot map
    % onto NNV's single-network linear reachability is NOT guessed at. Instead the
    % parser returns
    %     property.unsupported = true
    %     property.reason      = '<why>'
    % and (when possible) still fills a sound INPUT box so the caller can fall back
    % to falsification/PGD (a concrete counterexample is always a valid `sat`). The
    % runner treats `unsupported` as `unknown` for the reach verdict. Gated cases:
    % multiple networks (equal-to/isomorphic-to), declare-hidden, >1 input or output
    % tensor (multimodal), `!=`, nonlinear/arithmetic (`*`,`+`,`-` over variables),
    % partial indexing, and mixed input/output disjunctions.
    %
    % PERFORMANCE: the file is processed by STREAMING complete S-expression
    % statements (one (declare-network ...) / (assert ...) at a time) so the giant
    % image benchmarks (collins 96 MB, smart_turn 124 MB) never materialize a
    % whole-file token tree. Multi-network / multimodal files are gated from the
    % header alone -- we stop reading before the millions of input asserts.
    %
    % VNN-LIB 2.0 grammar handled (single network):
    %   (declare-network <name> (declare-input <Xname> <type> [shape])
    %                           (declare-output <Yname> <type> [shape]))
    %   (assert (<=|>=|<|>|== <Xname>[i,j,..] <const>))           input box
    %   (assert (and <X bounds> ...))                            input box (paired)
    %   (assert (or (and <X bounds>) (and <X bounds>) ...))      input OR of boxes
    %   (assert (<=|>=|<|> <Yname>[..] <const|Yname[..]>))        linear output
    %   (assert (or (and <Y lin> ...) <Y lin> ...))              disjunctive output
    %
    % Multi-dimensional indices are flattened ROW-MAJOR (ONNX/VNN-LIB order), the
    % same flat order load_vnnlib emits, so the runner's per-category `needReshape`
    % handling is unchanged. See VNNLIB2_SUPPORT_PLAN.md for the full design.

    fid = fopen(propertyFile, 'r');
    if fid < 0
        property = unsupported_property([], 'cannot open file');
        return;
    end
    closer = onCleanup(@() fclose(fid));

    st = init_state();
    depth = 0;              % running paren depth
    buf = '';               % current top-level statement being assembled
    started = false;        % has the current statement opened its first '(' ?

    while true
        line = fgetl(fid);
        if ~ischar(line), break; end
        line = strip_comment(line);
        if isempty(strtrim(line)), continue; end
        % scan char-by-char, emitting each balanced top-level (...) statement
        [st, depth, buf, started, done] = consume_chars(line, st, depth, buf, started);
        if done
            property = finalize(st);   % header/assert gated -> stop reading the rest
            return;
        end
    end

    property = finalize(st);
end

% =========================================================================
% Streaming statement assembler
% =========================================================================

function [st, depth, buf, started, done] = consume_chars(line, st, depth, buf, started)
    done = false;
    line = [line ' '];   % ensure a separator at line breaks
    for i = 1:numel(line)
        c = line(i);
        if c == '('
            depth = depth + 1; started = true;
        elseif c == ')'
            depth = depth - 1;
        end
        if started
            buf(end+1) = c; %#ok<AGROW>
        end
        if started && depth == 0
            stmt = strtrim(buf);
            buf = ''; started = false;
            if ~isempty(stmt)
                [st, done] = dispatch_statement(stmt, st);
                if done, return; end
            end
        end
    end
end

function [st, done] = dispatch_statement(stmt, st)
    done = false;
    head = statement_head(stmt);
    switch head
        case 'vnnlib-version'
            st.sawVersion = true;
        case 'declare-network'
            if st.sawNetwork
                st.unsupported = true;
                st.reason = 'multiple declare-network (multi-network 2.0 not supported)';
                done = true; return;
            end
            st.sawNetwork = true;
            st = parse_network_header(stmt, st);
            if st.unsupported
                done = true; return;   % gate from header alone
            end
        case 'assert'
            st = process_assert(stmt, st);
            if st.hardStop, done = true; return; end
        otherwise
            % ignore (set-logic, declare-const stragglers, comments-only, etc.)
    end
end

% =========================================================================
% Network header
% =========================================================================

function st = parse_network_header(stmt, st)
    node = parse_one_sexpr(stmt);
    inputs = {}; outputs = {}; hasEquiv = false; hasHidden = false;
    for k = 2:numel(node)
        c = node{k};
        if ~iscell(c) || isempty(c), continue; end
        h = char(c{1});
        switch h
            case 'declare-input',  inputs{end+1}  = c; %#ok<AGROW>
            case 'declare-output', outputs{end+1} = c; %#ok<AGROW>
            case 'declare-hidden', hasHidden = true;
            case {'equal-to','isomorphic-to'}, hasEquiv = true;
        end
    end
    if hasEquiv
        st.unsupported = true; st.reason = 'network equivalence (equal-to/isomorphic-to) -- multi-network not supported'; return;
    end
    if hasHidden
        st.unsupported = true; st.reason = 'declare-hidden intermediate constraints not supported'; return;
    end
    if numel(inputs) ~= 1 || numel(outputs) ~= 1
        st.unsupported = true;
        st.reason = sprintf('expected 1 input and 1 output tensor, found %d/%d (multimodal not supported)', numel(inputs), numel(outputs));
        return;
    end
    [st.inName, inShape]   = parse_io_decl(inputs{1});
    [st.outName, outShape] = parse_io_decl(outputs{1});
    st.inShape = inShape; st.outShape = outShape;
    st.inDim  = prod(inShape);
    st.outDim = prod(outShape);
    if isempty(st.inDim) || st.inDim < 1 || isempty(st.outDim) || st.outDim < 1
        st.unsupported = true; st.reason = 'could not parse input/output shape'; return;
    end
    st.lb = nan(st.inDim, 1, 'single');
    st.ub = nan(st.inDim, 1, 'single');
    st.haveBox = true;
end

% =========================================================================
% Assert processing
% =========================================================================

function st = process_assert(stmt, st)
    if ~st.sawNetwork || ~st.haveBox
        st.unsupported = true; st.reason = 'assert before declare-network'; st.hardStop = true; return;
    end

    % cheap, tree-free classification + unsupported-op scan on the raw statement
    refIn  = ref_present(stmt, st.inName);
    refOut = ref_present(stmt, st.outName);
    nonlin = has_unsupported_op_str(stmt);

    if refIn && refOut
        st.unsupported = true;
        st.reason = 'mixed input/output assertion (combined disjunction) not supported';
        st.hardStop = true; return;
    end

    if refIn
        if nonlin
            % nonlinear input region: cannot bound soundly -> gate (still record any
            % linear bounds we already have for PGD).
            st.unsupported = true;
            st.reason = 'nonlinear / arithmetic input constraint not soundly boundable';
            st.hardStop = true; return;
        end
        e = sexpr_body(parse_one_sexpr(stmt));   % the asserted expr
        head = node_head(e);
        if strcmp(head, 'or')
            [lbB, ubB, ok] = parse_input_or(e, st, st.lb, st.ub);
            if ~ok
                st.unsupported = true; st.reason = 'unsupported input OR form'; st.hardStop = true; return;
            end
            st.lbBoxes = lbB; st.ubBoxes = ubB; st.inputIsOr = true;
        else
            [st.lb, st.ub, ok] = apply_input_expr(e, st.lb, st.ub, st);
            if ~ok
                st.unsupported = true; st.reason = 'unsupported input constraint form'; st.hardStop = true; return;
            end
        end
    elseif refOut
        if nonlin
            st.unsupported = true;
            st.reason = 'nonlinear / arithmetic / != output property not soundly verifiable by linear reach';
            st.hardStop = true; return;
        end
        e = sexpr_body(parse_one_sexpr(stmt));
        [ast, ok] = build_output_assertion(e, st);
        if ~ok
            st.unsupported = true; st.reason = 'unsupported output property form'; st.hardStop = true; return;
        end
        [st.propHs, okacc] = accumulate_output(st.propHs, ast);
        if ~okacc
            st.unsupported = true;
            st.reason = 'output is a conjunction of multiple OR-clauses (product-of-sums) not representable as one disjunctive spec';
            st.hardStop = true; return;
        end
        st.sawOutput = true;
    else
        % references neither X nor Y -- ignore (constant tautology / set-info)
    end
end

% =========================================================================
% Finalization
% =========================================================================

function property = finalize(st)
    if ~st.sawVersion && ~st.sawNetwork
        property = unsupported_property([], 'no (vnnlib-version)/declare-network header -- not a 2.0 file'); return;
    end
    if st.unsupported
        box = [];
        if st.haveBox && ~st.inputIsOr && ~any(isnan(st.lb)) && ~any(isnan(st.ub))
            box = struct('lb', st.lb, 'ub', st.ub);
        end
        property = unsupported_property(box, st.reason);
        return;
    end
    if ~st.sawNetwork
        property = unsupported_property([], 'no declare-network'); return;
    end

    if st.inputIsOr
        for i = 1:numel(st.lbBoxes)
            if any(isnan(st.lbBoxes{i})) || any(isnan(st.ubBoxes{i}))
                property = unsupported_property([], 'input OR box has unbounded dimension'); return;
            end
        end
        property.lb = st.lbBoxes;
        property.ub = st.ubBoxes;
    else
        if any(isnan(st.lb)) || any(isnan(st.ub))
            property = unsupported_property([], 'input box has an unbounded dimension (under-constrained)'); return;
        end
        property.lb = st.lb;
        property.ub = st.ub;
    end

    if ~st.sawOutput || isempty(st.propHs)
        property.unsupported = true;
        property.reason = 'no output property assertion found';
        property.prop = {};
        return;
    end
    property.prop = st.propHs;
    property.unsupported = false;
end

function st = init_state()
    st = struct();
    st.sawVersion = false;
    st.sawNetwork = false;
    st.unsupported = false;
    st.reason = '';
    st.hardStop = false;
    st.inName = ''; st.outName = '';
    st.inShape = []; st.outShape = [];
    st.inDim = []; st.outDim = [];
    st.lb = []; st.ub = []; st.haveBox = false;
    st.lbBoxes = {}; st.ubBoxes = {}; st.inputIsOr = false;
    st.propHs = {}; st.sawOutput = false;
end

% =========================================================================
% Cheap raw-string helpers (no tree build)
% =========================================================================

function s = strip_comment(line)
    ci = strfind(line, ';');
    if isempty(ci), s = line; else, s = line(1:ci(1)-1); end
end

function h = statement_head(stmt)
    % first token after the opening '('
    h = '';
    t = regexp(stmt, '^\(\s*([A-Za-z0-9_\-!]+)', 'tokens', 'once');
    if ~isempty(t), h = t{1}; end
end

function tf = ref_present(stmt, name)
    % does the variable `name` appear as a token (name[...] or bare name)? short-circuits.
    if isempty(name), tf = false; return; end
    pat = ['(?<![A-Za-z0-9_])' regexptranslate('escape', name) '(?![A-Za-z0-9_])'];
    tf = ~isempty(regexp(stmt, pat, 'once'));
end

function tf = has_unsupported_op_str(stmt)
    % arithmetic operators applied to variables, or `!=`. Constants like -0.5 are
    % NOT `(- ...)` operators, so we look for the operator-in-call form `(op `.
    tf = ~isempty(regexp(stmt, '\(\s*[\*\+/]', 'once')) || ...
         ~isempty(regexp(stmt, '\(\s*-\s', 'once')) || ...
         contains(stmt, '!=');
end

% =========================================================================
% Tiny S-expression parser (per-statement; operates on ONE statement string)
% =========================================================================

function node = parse_one_sexpr(stmt)
    s = squeeze_bracket_space(stmt);
    toks = tokenize(s);
    i = 1;
    if i <= numel(toks) && strcmp(toks{i}, '(')
        node = parse_node(toks, i);
    else
        node = {};
    end
end

function body = sexpr_body(node)
    % the asserted expression: (assert <expr>) -> <expr>
    if iscell(node) && numel(node) >= 2
        body = node{2};
    else
        body = {};
    end
end

function s = squeeze_bracket_space(s)
    % remove whitespace inside [...] so "X[0, 0, 0, 3]" / "[1, 1, 1, 5]" tokenize whole
    out = char(zeros(1, numel(s)));
    n = 0; depth = 0;
    for i = 1:numel(s)
        c = s(i);
        if c == '[', depth = depth + 1; elseif c == ']', depth = max(0, depth - 1); end
        if depth > 0 && (c == ' ' || c == sprintf('\t')), continue; end
        n = n + 1; out(n) = c;
    end
    s = out(1:n);
end

function toks = tokenize(s)
    s = strrep(s, '(', ' ( ');
    s = strrep(s, ')', ' ) ');
    toks = regexp(strtrim(s), '\s+', 'split');
    toks = toks(~cellfun(@isempty, toks));
end

function node = parse_node(toks, i)
    % toks{i} == '('
    node = {};
    i = i + 1;
    while i <= numel(toks)
        t = toks{i};
        if strcmp(t, '(')
            [child, i] = parse_node_r(toks, i);
            node{end+1} = child; %#ok<AGROW>
        elseif strcmp(t, ')')
            return;
        else
            node{end+1} = t; %#ok<AGROW>
            i = i + 1;
        end
    end
end

function [node, i] = parse_node_r(toks, i)
    node = {};
    i = i + 1;
    while i <= numel(toks)
        t = toks{i};
        if strcmp(t, '(')
            [child, i] = parse_node_r(toks, i);
            node{end+1} = child; %#ok<AGROW>
        elseif strcmp(t, ')')
            i = i + 1;
            return;
        else
            node{end+1} = t; %#ok<AGROW>
            i = i + 1;
        end
    end
end

function h = node_head(e)
    h = '';
    if iscell(e) && ~isempty(e) && (ischar(e{1}) || isstring(e{1}))
        h = char(e{1});
    end
end

% =========================================================================
% Declaration parsing
% =========================================================================

function [name, shape] = parse_io_decl(decl)
    % decl = {'declare-input'/'declare-output', name, type, '[d0,d1,..]'}
    name = ''; shape = [];
    if numel(decl) < 2, return; end
    name = char(decl{2});
    for k = numel(decl):-1:3
        if ischar(decl{k}) || isstring(decl{k})
            tk = char(decl{k});
            if startsWith(tk, '[') && endsWith(tk, ']')
                inner = tk(2:end-1);
                if isempty(inner)
                    shape = 1;            % scalar
                else
                    shape = str2double(regexp(inner, ',', 'split'));
                end
                return;
            end
        end
    end
end

% =========================================================================
% Index flattening
% =========================================================================

function nm = var_basename(tok)
    nm = ''; tok = char(tok);
    b = strfind(tok, '[');
    if ~isempty(b), cand = tok(1:b(1)-1); else, cand = tok; end
    if ~isempty(cand) && (isletter(cand(1)) || cand(1) == '_'), nm = cand; end
end

function [flat, ok] = flatten_ref(tok, name, shape)
    % 'X[0,0,0,3]' with shape [1 1 1 5] -> flat 0-based row-major index (here 3).
    flat = []; ok = false; tok = char(tok);
    b = strfind(tok, '[');
    if isempty(b)
        if prod(shape) == 1, flat = 0; ok = true; end
        return;
    end
    if ~strcmp(tok(1:b(1)-1), name), return; end
    e = strfind(tok, ']');
    if isempty(e), return; end
    inner = tok(b(1)+1:e(end)-1);
    idx = str2double(regexp(inner, ',', 'split'));
    if any(isnan(idx)), return; end
    if numel(idx) ~= numel(shape), return; end   % partial indexing illegal in 2.0
    if any(idx < 0) || any(idx >= shape), return; end
    f = 0;
    for k = 1:numel(shape)
        f = f * shape(k) + idx(k);
    end
    flat = f; ok = true;
end

% =========================================================================
% Input constraints
% =========================================================================

function [lb, ub, ok] = apply_input_expr(e, lb, ub, st)
    ok = true;
    if ~iscell(e) || isempty(e), ok = false; return; end
    head = char(e{1});
    if strcmp(head, 'and')
        for k = 2:numel(e)
            [lb, ub, ok] = apply_input_expr(e{k}, lb, ub, st);
            if ~ok, return; end
        end
        return;
    end
    if ~ismember(head, {'<=','>=','<','>','=='}), ok = false; return; end
    if numel(e) ~= 3, ok = false; return; end
    var = char(e{2}); val = e{3};
    if ~ischar(val) && ~isstring(val), ok = false; return; end
    if ~isempty(var_basename(val)), ok = false; return; end   % var-var input not a box
    c = single(str2double(char(val)));
    if isnan(c), ok = false; return; end
    [flat, okf] = flatten_ref(var, st.inName, st.inShape);
    if ~okf, ok = false; return; end
    idx = flat + 1;
    % Strict (<,>) are treated as non-strict (<=,>=), matching the battle-tested 1.0
    % parser (load_vnnlib.m process_constraint/process_input_constraint). This is sound
    % for the UNSAT proof verdict (the box/unsafe region is over-approximated). For SAT
    % the only divergence is a witness EXACTLY on the boundary, which (a) is guarded by
    % the onnxruntime witness replay and (b) falls within VNN-COMP's 1e-4 constraint
    % tolerance, so the official checker treats it identically.
    switch head
        case {'>=','>'}, lb(idx) = c;
        case {'<=','<'}, ub(idx) = c;
        case '==',       lb(idx) = c; ub(idx) = c;
    end
end

function [lbBoxes, ubBoxes, ok] = parse_input_or(e, st, lb0, ub0)
    lbBoxes = {}; ubBoxes = {}; ok = true;
    for k = 2:numel(e)
        d = e{k};
        lb = lb0; ub = ub0;
        [lb, ub, okd] = apply_input_expr(d, lb, ub, st);
        if ~okd, ok = false; return; end
        lbBoxes{end+1} = lb; %#ok<AGROW>
        ubBoxes{end+1} = ub; %#ok<AGROW>
    end
    if isempty(lbBoxes), ok = false; end
end

% =========================================================================
% Output property (unsafe region) -> HalfSpace
% =========================================================================

function [ast, ok] = build_output_assertion(e, st)
    ast = struct('Hg', HalfSpace.empty); ok = true;
    if ~iscell(e) || isempty(e), ok = false; return; end
    head = char(e{1});
    if strcmp(head, 'or')
        Hs = HalfSpace.empty;
        for k = 2:numel(e)
            [sub, oks] = output_conjunction(e{k}, st);
            if ~oks, ok = false; return; end
            Hs(end+1,1) = sub; %#ok<AGROW>
        end
        ast.Hg = Hs;
    else
        [Hs, oks] = output_conjunction(e, st);
        if ~oks, ok = false; return; end
        ast.Hg = Hs;
    end
end

function [hs, ok] = output_conjunction(e, st)
    ok = true; G = []; g = [];
    head = node_head(e);
    if strcmp(head, 'and')
        for k = 2:numel(e)
            [Gr, gr, okr] = output_rows(e{k}, st);
            if ~okr, ok = false; hs = HalfSpace.empty; return; end
            G = [G; Gr]; g = [g; gr]; %#ok<AGROW>
        end
    else
        [G, g, okr] = output_rows(e, st);
        if ~okr, ok = false; hs = HalfSpace.empty; return; end
    end
    hs = HalfSpace(G, g);
end

function [G, g, ok] = output_rows(e, st)
    ok = true; G = []; g = [];
    if ~iscell(e) || numel(e) ~= 3, ok = false; return; end
    head = char(e{1});
    if ~ismember(head, {'<=','>=','<','>','=='}), ok = false; return; end
    [row, gval, okr] = linear_two_terms(head, e{2}, e{3}, st);
    if ~okr, ok = false; return; end
    if strcmp(head, '==')
        G = [row; -row]; g = [gval; -gval];
    else
        G = row; g = gval;
    end
end

function [row, gval, ok] = linear_two_terms(op, lhs, rhs, st)
    ok = true; row = zeros(1, st.outDim); gval = 0;
    [li, lc, okl] = term_value(lhs, st);
    [ri, rc, okr] = term_value(rhs, st);
    if ~okl || ~okr, ok = false; return; end
    coef = zeros(1, st.outDim);
    cst = 0;
    if ~isempty(li), coef(li) = coef(li) + 1; else, cst = cst - lc; end
    if ~isempty(ri), coef(ri) = coef(ri) - 1; else, cst = cst + rc; end
    if any(strcmp(op, {'<=','<'}))
        row = coef; gval = single(cst);
    else
        row = -coef; gval = single(-cst);
    end
end

function [idx, c, ok] = term_value(t, st)
    idx = []; c = 0; ok = true; t = char(t);
    nm = var_basename(t);
    if isempty(nm)
        c = single(str2double(t));
        if isnan(c), ok = false; end
        return;
    end
    if ~strcmp(nm, st.outName), ok = false; return; end
    [flat, okf] = flatten_ref(t, st.outName, st.outShape);
    if ~okf, ok = false; return; end
    idx = flat + 1;
end

function [propHs, ok] = accumulate_output(propHs, ast)
    % Combine successive output assertions into ONE prop cell. The runner verifies a
    % single-input property against prop{1}.Hg only (run_vnncomp_instance.m: the
    % length(prop)==1 branches), so we must NEVER leave multiple cells -- a trailing
    % conjunct after an (or ...) has to be AND-ed INTO every disjunct, not appended as
    % a second cell that the runner would silently ignore (-150 false verdict). This
    % mirrors the hardened 1.0 parser's distribution (load_vnnlib.m ~116-153).
    ok = true;
    if isempty(propHs), propHs = {ast}; return; end
    last = propHs{end};
    nL = numel(last.Hg); nA = numel(ast.Hg);
    if nL == 1 && nA == 1
        % both scalar conjuncts -> AND-merge rows into one polytope
        last.Hg.G = [last.Hg.G; ast.Hg.G];
        last.Hg.g = [last.Hg.g; ast.Hg.g];
        propHs{end} = last;
    elseif nL > 1 && nA == 1
        % running region is an OR; new standalone conjunct ANDs into EVERY disjunct
        for di = 1:nL
            last.Hg(di).G = [last.Hg(di).G; ast.Hg.G];
            last.Hg(di).g = [last.Hg(di).g; ast.Hg.g];
        end
        propHs{end} = last;
    elseif nL == 1 && nA > 1
        % running region is a scalar conjunct asserted BEFORE this OR; AND the
        % conjunct's rows into each disjunct of the new OR and keep the OR
        for di = 1:nA
            ast.Hg(di).G = [ast.Hg(di).G; last.Hg.G];
            ast.Hg(di).g = [ast.Hg(di).g; last.Hg.g];
        end
        propHs{end} = ast;
    else
        % OR conjoined with OR: a product-of-sums NNV's single-disjunctive-spec
        % contract cannot represent (the 1.0 parser errors here) -> gate (sound).
        ok = false;
    end
end

% =========================================================================
% Unsupported-property helpers
% =========================================================================

function property = unsupported_property(box, reason)
    property = struct();
    property.unsupported = true;
    property.reason = reason;
    property.prop = {};
    if ~isempty(box) && isfield(box, 'lb')
        property.lb = box.lb; property.ub = box.ub;
    else
        property.lb = []; property.ub = [];
    end
end
