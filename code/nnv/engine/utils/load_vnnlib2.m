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
    % and (when possible) still records a sound INPUT box in property.lb/ub. NOTE: the
    % current run_vnncomp_instance.m emits `unknown` for ANY unsupported property; the
    % recorded box is reserved for a FUTURE falsification/PGD fallback (plan Phase 3b)
    % that could still find a concrete `sat` for a gated-but-linear-input case. Gated:
    % three-or-more networks, isomorphic-to, declare-hidden, >1 input or output
    % tensor (multimodal), `!=`, nonlinear/arithmetic (`*`,`+`,`-` over variables),
    % partial indexing, and mixed input/output disjunctions.
    %
    % PHASE 2 (multi-network, `equal-to` ONLY -- plan section 3.1): when the file
    % declares exactly TWO networks and the second is `(equal-to <first>)` (the
    % monotonic_acasxu_2026 shape: the SAME onnx evaluated on two coupled inputs),
    % the file is PARSED instead of gated. The result keeps the single-network
    % contract fields EMPTY (property.lb/ub = [], property.prop = {}) so legacy
    % consumers cannot misuse them, and adds
    %     property.multinet.names     -> {name_f, name_g}            (cell 1x2)
    %     property.multinet.inShapes  -> {inShape_f, inShape_g}      (must match)
    %     property.multinet.outShapes -> {outShape_f, outShape_g}    (must match)
    %     property.multinet.equivKind -> 'equal'
    %     property.multinet.jointLb   -> box over the STACKED input [X_f; X_g] (2n x 1)
    %     property.multinet.jointUb      (g-block bounds that only follow from the
    %                                     coupling rows are closed by sound interval
    %                                     propagation; see propagate_joint_bounds)
    %     property.multinet.jointC    -> cross/within-network input coupling rows
    %     property.multinet.jointd       jointC * x <= jointd over the CONCRETE
    %                                     stacked input x = [X_f; X_g] (NOT over star
    %                                     predicates -- verify_multinet maps x to the
    %                                     Star predicate alpha; `==` couplings are
    %                                     emitted as two-row inequality pairs)
    %     property.multinet.crossProp -> HalfSpace array over the STACKED output
    %                                     [Y_f; Y_g]: G*y <= g is the UNSAFE region,
    %                                     same polarity convention as .prop{n}.Hg
    % property.unsupported stays true (the default) unless EVERY assert parsed
    % cleanly into the joint box / jointC / crossProp AND the joint box closed.
    % Verification of this structure lives in engine/utils/verify_multinet.m.
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
                % Phase 2: a SECOND declare-network is parsed (not gated) iff it
                % is `(equal-to <first>)` with matching shapes; anything else
                % (third network, isomorphic-to, multimodal, ...) stays gated.
                st = parse_second_network(stmt, st);
                if st.unsupported
                    done = true; return;
                end
            else
                st.sawNetwork = true;
                st = parse_network_header(stmt, st);
                if st.unsupported
                    done = true; return;   % gate from header alone
                end
            end
        case 'assert'
            st.sawAssert = true;
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
    if numel(node) >= 2 && (ischar(node{2}) || isstring(node{2}))
        st.netName = char(node{2});   % needed to validate `(equal-to <name>)` (Phase 2)
    end
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

function st = parse_second_network(stmt, st)
    % Phase 2 (multi-network `equal-to` ONLY, plan section 3.1): accept a SECOND
    % declare-network iff it is `(equal-to <first-network-name>)` with exactly one
    % input and one output tensor whose shapes MATCH the first network (the 2.0
    % standard requires matching element types and shapes for equal-to). On
    % success the parser switches to multinet mode: subsequent asserts fill the
    % joint stacked box/coupling rows/cross output property instead of the
    % single-network fields. EVERY other shape of file stays gated (sound).
    if st.multinet
        st.unsupported = true;
        st.reason = 'more than two declare-network (only an equal-to pair is supported)';
        return;
    end
    if st.sawAssert
        % The 2.0 grammar puts all asserts AFTER the network declarations. An
        % assert seen BEFORE the second declare-network was processed (or
        % silently ignored, if it referenced the then-undeclared network) in
        % single-network mode, so continuing could DROP constraints -- a sat
        % witness ignoring a dropped input constraint fails the official
        % witness replay (-150). Gate instead.
        st.unsupported = true;
        st.reason = 'assert appears before the second declare-network';
        return;
    end
    node = parse_one_sexpr(stmt);
    name2 = '';
    if numel(node) >= 2 && (ischar(node{2}) || isstring(node{2}))
        name2 = char(node{2});
    end
    inputs = {}; outputs = {}; equivs = {}; hasHidden = false;
    for k = 2:numel(node)
        c = node{k};
        if ~iscell(c) || isempty(c), continue; end
        h = node_head(c);
        switch h
            case 'declare-input',  inputs{end+1}  = c; %#ok<AGROW>
            case 'declare-output', outputs{end+1} = c; %#ok<AGROW>
            case 'declare-hidden', hasHidden = true;
            case {'equal-to','isomorphic-to'}, equivs{end+1} = c; %#ok<AGROW>
        end
    end
    if hasHidden
        st.unsupported = true; st.reason = 'declare-hidden intermediate constraints not supported'; return;
    end
    if numel(equivs) ~= 1
        st.unsupported = true; st.reason = 'second declare-network must carry exactly one equivalence relation (equal-to)'; return;
    end
    eq = equivs{1};
    if ~strcmp(node_head(eq), 'equal-to')
        st.unsupported = true; st.reason = 'isomorphic-to (different weights) not supported -- only equal-to'; return;
    end
    if numel(eq) < 2 || ~(ischar(eq{2}) || isstring(eq{2})) || isempty(st.netName) || ~strcmp(char(eq{2}), st.netName)
        st.unsupported = true; st.reason = 'equal-to must reference the first declared network'; return;
    end
    if numel(inputs) ~= 1 || numel(outputs) ~= 1
        st.unsupported = true;
        st.reason = sprintf('expected 1 input and 1 output tensor per network, found %d/%d (multimodal not supported)', numel(inputs), numel(outputs));
        return;
    end
    [st.inName2, inShape2]   = parse_io_decl(inputs{1});
    [st.outName2, outShape2] = parse_io_decl(outputs{1});
    st.inShape2 = inShape2; st.outShape2 = outShape2;
    if ~isequal(inShape2, st.inShape) || ~isequal(outShape2, st.outShape)
        st.unsupported = true; st.reason = 'equal-to networks must have identical input/output shapes'; return;
    end
    % the four tensor names (and the two net names) must be distinct, otherwise
    % the raw-string assert classification (ref_present) is ambiguous
    nm = {st.inName, st.outName, st.inName2, st.outName2};
    if any(cellfun(@isempty, nm)) || numel(unique(nm)) ~= 4 || isempty(name2) || strcmp(name2, st.netName)
        st.unsupported = true; st.reason = 'ambiguous / duplicate tensor names across the two networks'; return;
    end
    % switch to multinet mode: stacked order is [X_f; X_g] / [Y_f; Y_g]
    st.multinet = true;
    st.netName2 = name2;
    st.jointLb = nan(2*st.inDim, 1, 'single');
    st.jointUb = nan(2*st.inDim, 1, 'single');
    st.jointC = zeros(0, 2*st.inDim, 'single');
    st.jointd = zeros(0, 1, 'single');
    st.crossCells = {};
end

% =========================================================================
% Assert processing
% =========================================================================

function st = process_assert(stmt, st)
    if st.multinet
        st = process_assert_multinet(stmt, st);
        return;
    end
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
% Assert processing -- multi-network (equal-to) mode (Phase 2)
% =========================================================================

function st = process_assert_multinet(stmt, st)
    % Multinet twin of process_assert. Stacked variable order is fixed as
    % [X_f; X_g] for inputs and [Y_f; Y_g] for outputs (f = first declared
    % network). Any form that does not map EXACTLY onto
    %   - a var-const bound on the joint box,
    %   - a two-term linear coupling row jointC*x <= jointd, or
    %   - a linear (possibly disjunctive) HalfSpace over [Y_f; Y_g]
    % gates the whole property (unsupported -> runner says `unknown`).
    refIn  = ref_present(stmt, st.inName)  || ref_present(stmt, st.inName2);
    refOut = ref_present(stmt, st.outName) || ref_present(stmt, st.outName2);
    nonlin = has_unsupported_op_str(stmt);

    if refIn && refOut
        st.unsupported = true;
        st.reason = 'mixed input/output assertion in multi-network property not supported';
        st.hardStop = true; return;
    end

    if refIn
        if nonlin
            st.unsupported = true;
            st.reason = 'nonlinear / arithmetic multi-network input constraint not soundly representable';
            st.hardStop = true; return;
        end
        e = sexpr_body(parse_one_sexpr(stmt));
        [st, ok] = apply_mn_input_expr(e, st);
        if ~ok
            st.unsupported = true; st.reason = 'unsupported multi-network input constraint form'; st.hardStop = true; return;
        end
    elseif refOut
        if nonlin
            st.unsupported = true;
            st.reason = 'nonlinear / arithmetic / != multi-network output property not soundly verifiable by linear reach';
            st.hardStop = true; return;
        end
        e = sexpr_body(parse_one_sexpr(stmt));
        [ast, ok] = mn_build_output_assertion(e, st);
        if ~ok
            st.unsupported = true; st.reason = 'unsupported multi-network output property form'; st.hardStop = true; return;
        end
        [st.crossCells, okacc] = accumulate_output(st.crossCells, ast);
        if ~okacc
            st.unsupported = true;
            st.reason = 'multi-network output is a product-of-sums not representable as one disjunctive spec';
            st.hardStop = true; return;
        end
        st.sawCrossOutput = true;
    else
        % references neither network's tensors -- ignore (constant tautology)
    end
end

function [st, ok] = apply_mn_input_expr(e, st)
    % Multinet twin of apply_input_expr. Handles
    %   (and ...)                          recursion
    %   (op  <inVar> <const>)              joint box bound (strict -> non-strict,
    %                                      same soundness note as apply_input_expr)
    %   (op  <inVar> <inVar>)              linear coupling row(s) over the CONCRETE
    %                                      stacked input x = [X_f; X_g]:
    %       (<= A B), (< A B)   ->  x(A) - x(B) <= 0      (one row:  +1@A, -1@B)
    %       (>= A B), (> A B)   -> -x(A) + x(B) <= 0      (one row:  -1@A, +1@B)
    %       (== A B)            ->  both rows above        (two-row inequality pair)
    % NOTE: jointC/jointd constrain x DIRECTLY (concrete input space), NOT the
    % Star predicate alpha. verify_multinet performs the x -> alpha substitution
    % (x = c + G*alpha) when building the joint input Star; the falsifier checks
    % jointC*x <= jointd on concrete samples. Parser and verifier MUST stay
    % consistent on this convention.
    ok = true;
    if ~iscell(e) || isempty(e), ok = false; return; end
    head = char(e{1});
    if strcmp(head, 'and')
        for k = 2:numel(e)
            [st, ok] = apply_mn_input_expr(e{k}, st);
            if ~ok, return; end
        end
        return;
    end
    if ~ismember(head, {'<=','>=','<','>','=='}), ok = false; return; end
    if numel(e) ~= 3, ok = false; return; end
    lhs = e{2}; rhs = e{3};
    if ~(ischar(lhs) || isstring(lhs)) || ~(ischar(rhs) || isstring(rhs))
        ok = false; return;   % nested expression -- not a plain bound/coupling
    end
    [li, lok] = mn_input_index(lhs, st);
    if ~lok, ok = false; return; end   % lhs must be an input ref (mirrors 1-net path)
    if isempty(var_basename(rhs))
        % var-const -> joint box (strict <,> treated as <=,>= like apply_input_expr)
        c = single(str2double(char(rhs)));
        if isnan(c), ok = false; return; end
        switch head
            case {'>=','>'}, st.jointLb(li) = c;
            case {'<=','<'}, st.jointUb(li) = c;
            case '==',       st.jointLb(li) = c; st.jointUb(li) = c;
        end
    else
        % var-var -> coupling row(s); both sides must resolve to input dims
        [ri, rok] = mn_input_index(rhs, st);
        if ~rok, ok = false; return; end
        row = zeros(1, 2*st.inDim, 'single');
        switch head
            case {'<=','<'}
                row(li) = row(li) + 1; row(ri) = row(ri) - 1;       % x(li) - x(ri) <= 0
                st.jointC = [st.jointC; row];
                st.jointd = [st.jointd; single(0)];
            case {'>=','>'}
                row(li) = row(li) - 1; row(ri) = row(ri) + 1;       % x(ri) - x(li) <= 0
                st.jointC = [st.jointC; row];
                st.jointd = [st.jointd; single(0)];
            case '=='
                row(li) = row(li) + 1; row(ri) = row(ri) - 1;
                st.jointC = [st.jointC; row; -row];                 % == as two-row pair
                st.jointd = [st.jointd; single(0); single(0)];
        end
    end
end

function [idx, ok] = mn_input_index(tok, st)
    % resolve an input variable reference to its 1-based STACKED index in
    % x = [X_f; X_g] (f-block first, then the g-block offset by inDim)
    idx = []; ok = false; tok = char(tok);
    nm = var_basename(tok);
    if isempty(nm), return; end
    if strcmp(nm, st.inName)
        [flat, okf] = flatten_ref(tok, st.inName, st.inShape);
        if ~okf, return; end
        idx = flat + 1; ok = true;
    elseif strcmp(nm, st.inName2)
        [flat, okf] = flatten_ref(tok, st.inName2, st.inShape2);
        if ~okf, return; end
        idx = st.inDim + flat + 1; ok = true;
    end
end

function [ast, ok] = mn_build_output_assertion(e, st)
    % multinet twin of build_output_assertion (same or/and shapes)
    ast = struct('Hg', HalfSpace.empty); ok = true;
    if ~iscell(e) || isempty(e), ok = false; return; end
    head = char(e{1});
    if strcmp(head, 'or')
        Hs = HalfSpace.empty;
        for k = 2:numel(e)
            [sub, oks] = mn_output_conjunction(e{k}, st);
            if ~oks, ok = false; return; end
            Hs(end+1,1) = sub; %#ok<AGROW>
        end
        ast.Hg = Hs;
    else
        [Hs, oks] = mn_output_conjunction(e, st);
        if ~oks, ok = false; return; end
        ast.Hg = Hs;
    end
end

function [hs, ok] = mn_output_conjunction(e, st)
    % multinet twin of output_conjunction
    ok = true; G = []; g = [];
    head = node_head(e);
    if strcmp(head, 'and')
        for k = 2:numel(e)
            [Gr, gr, okr] = mn_output_rows(e{k}, st);
            if ~okr, ok = false; hs = HalfSpace.empty; return; end
            G = [G; Gr]; g = [g; gr]; %#ok<AGROW>
        end
    else
        [G, g, okr] = mn_output_rows(e, st);
        if ~okr, ok = false; hs = HalfSpace.empty; return; end
    end
    hs = HalfSpace(G, g);
end

function [G, g, ok] = mn_output_rows(e, st)
    % multinet twin of output_rows (== expands to a two-row pair)
    ok = true; G = []; g = [];
    if ~iscell(e) || numel(e) ~= 3, ok = false; return; end
    head = char(e{1});
    if ~ismember(head, {'<=','>=','<','>','=='}), ok = false; return; end
    [row, gval, okr] = mn_linear_two_terms(head, e{2}, e{3}, st);
    if ~okr, ok = false; return; end
    if strcmp(head, '==')
        G = [row; -row]; g = [gval; -gval];
    else
        G = row; g = gval;
    end
end

function [row, gval, ok] = mn_linear_two_terms(op, lhs, rhs, st)
    % Multinet twin of linear_two_terms over the STACKED output y = [Y_f; Y_g]
    % (width 2*outDim). POLARITY DERIVATION (the -150-critical part), mirroring
    % linear_two_terms exactly:
    %
    %   VNN-LIB asserts encode the SAT/counterexample (= unsafe) region ITSELF:
    %   a `sat` witness must make every assert TRUE, and the official checker
    %   REPLAYS the asserts on the witness. So the asserted comparison maps
    %   DIRECTLY into HalfSpace G*y <= g -- no negation. (Negating here, i.e.
    %   encoding the complement, would flip sat/unsat: a "witness" would be
    %   rejected by the replay and a -150 wrong verdict would follow.)
    %
    %   (op lhs rhs) is first normalized to  coef*y ? cst  with lhs terms +1,
    %   rhs terms -1 (constants folded into cst with opposite signs):
    %     op in {<=,<}:  lhs - rhs <= 0  ->  row =  coef, g =  cst
    %     op in {>=,>}:  lhs - rhs >= 0  ->  row = -coef, g = -cst
    %   Strict <,> are relaxed to <=,>= (superset of the true unsafe region:
    %   sound for the UNSAT proof; the boundary-witness caveat is the same as
    %   in the single-network parser and is guarded by strict-interior checks
    %   in verify_multinet's falsifier).
    %
    %   Worked example (monotonic_acasxu): assert (< Y_f[3] Y_g[3]) ->
    %   unsafe = { Y_f[3] - Y_g[3] <= 0 }, i.e. row has +1 at the Y_f[3] column
    %   (stacked col 4) and -1 at the Y_g[3] column (stacked col outDim+4),
    %   g = 0. An `unsat` verdict then proves NO coupled input pair reaches
    %   Y_f[3] <= Y_g[3], hence none reaches the strict Y_f[3] < Y_g[3] either
    %   (the asserted region is unreachable -> property proved).
    ok = true; row = zeros(1, 2*st.outDim); gval = 0;
    [li, lc, okl] = mn_term_value(lhs, st);
    [ri, rc, okr] = mn_term_value(rhs, st);
    if ~okl || ~okr, ok = false; return; end
    coef = zeros(1, 2*st.outDim);
    cst = 0;
    if ~isempty(li), coef(li) = coef(li) + 1; else, cst = cst - lc; end
    if ~isempty(ri), coef(ri) = coef(ri) - 1; else, cst = cst + rc; end
    if any(strcmp(op, {'<=','<'}))
        row = coef; gval = single(cst);
    else
        row = -coef; gval = single(-cst);
    end
end

function [idx, c, ok] = mn_term_value(t, st)
    % multinet twin of term_value: resolve an output ref into the STACKED
    % [Y_f; Y_g] order (g-block offset by outDim), or a numeric constant
    idx = []; c = 0; ok = true; t = char(t);
    nm = var_basename(t);
    if isempty(nm)
        c = single(str2double(t));
        if isnan(c), ok = false; end
        return;
    end
    if strcmp(nm, st.outName)
        [flat, okf] = flatten_ref(t, st.outName, st.outShape);
        if ~okf, ok = false; return; end
        idx = flat + 1;
    elseif strcmp(nm, st.outName2)
        [flat, okf] = flatten_ref(t, st.outName2, st.outShape2);
        if ~okf, ok = false; return; end
        idx = st.outDim + flat + 1;
    else
        ok = false;
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
        % never record a single-net box for a multinet file: it would describe
        % only the first network's input and could mislead a future PGD fallback
        if ~st.multinet && st.haveBox && ~st.inputIsOr && ~any(isnan(st.lb)) && ~any(isnan(st.ub))
            box = struct('lb', st.lb, 'ub', st.ub);
        end
        property = unsupported_property(box, st.reason);
        return;
    end
    if ~st.sawNetwork
        property = unsupported_property([], 'no declare-network'); return;
    end
    if st.multinet
        property = finalize_multinet(st);
        return;
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
    st.sawAssert = false;
    st.unsupported = false;
    st.reason = '';
    st.hardStop = false;
    st.netName = '';
    st.inName = ''; st.outName = '';
    st.inShape = []; st.outShape = [];
    st.inDim = []; st.outDim = [];
    st.lb = []; st.ub = []; st.haveBox = false;
    st.lbBoxes = {}; st.ubBoxes = {}; st.inputIsOr = false;
    st.propHs = {}; st.sawOutput = false;
    % --- Phase 2: multi-network (equal-to) state ---
    st.multinet = false;
    st.netName2 = '';
    st.inName2 = ''; st.outName2 = '';
    st.inShape2 = []; st.outShape2 = [];
    st.jointLb = []; st.jointUb = [];
    st.jointC = []; st.jointd = [];
    st.crossCells = {}; st.sawCrossOutput = false;
end

function property = finalize_multinet(st)
    % Phase 2 finalization. UNSUPPORTED-BY-DEFAULT discipline: every gating
    % assert already set st.unsupported (handled by the caller before we get
    % here), so at this point all asserts parsed cleanly; what remains is to
    % (a) require a cross-network output property,
    % (b) CLOSE the joint box -- the g-block is typically bounded only THROUGH
    %     the coupling rows (X_g == X_f componentwise, X_g[0] <= X_f[0], ...),
    %     so propagate interval bounds through jointC*x <= jointd, and
    % (c) refuse anything still unbounded / infeasible.
    if ~st.sawCrossOutput || isempty(st.crossCells)
        property = unsupported_property([], 'multi-network: no output property assertion found');
        return;
    end
    if numel(st.crossCells) ~= 1
        % accumulate_output keeps exactly one cell on success; defensive only
        property = unsupported_property([], 'multi-network: output property is not a single disjunctive spec');
        return;
    end
    [jlb, jub] = propagate_joint_bounds(st.jointLb, st.jointUb, st.jointC, st.jointd);
    if any(isnan(jlb)) || any(isnan(jub))
        property = unsupported_property([], 'multi-network: joint input box has an unbounded dimension (under-constrained)');
        return;
    end
    if any(jlb > jub)
        % constraints are infeasible; technically vacuously unsat, but we gate
        % conservatively instead of emitting a verdict from an empty region
        property = unsupported_property([], 'multi-network: joint input constraints are infeasible (empty box)');
        return;
    end
    property = struct();
    property.lb = [];     % single-network contract fields left EMPTY on purpose:
    property.ub = [];     % legacy single-net consumers must not see a usable box
    property.prop = {};
    mn = struct();
    mn.names = {st.netName, st.netName2};
    mn.inShapes = {st.inShape, st.inShape2};
    mn.outShapes = {st.outShape, st.outShape2};
    mn.equivKind = 'equal';
    mn.jointLb = jlb;
    mn.jointUb = jub;
    mn.jointC = st.jointC;
    mn.jointd = st.jointd;
    mn.crossProp = st.crossCells{1}.Hg;
    property.multinet = mn;
    property.unsupported = false;
    property.reason = '';
end

function [lb, ub] = propagate_joint_bounds(lb, ub, C, d)
    % Sound interval tightening of the joint box through the two-term coupling
    % rows of C*x <= d (the parser only emits rows with exactly two +/-1
    % entries; the algebra below is generic anyway). For a row
    %       a_i*x_i + a_j*x_j <= d_r
    % every feasible point satisfies a_i*x_i <= d_r - a_j*x_j, hence
    %       a_i*x_i <= d_r - min_{x_j in [lb_j,ub_j]}(a_j*x_j)
    % with min(a_j*x_j) = a_j*lb_j if a_j>0, else a_j*ub_j. Dividing by a_i:
    %       a_i > 0  ->  x_i <= (d_r - min(a_j*x_j))/a_i   (an UPPER bound)
    %       a_i < 0  ->  x_i >= (d_r - min(a_j*x_j))/a_i   (a  LOWER bound)
    % NaN (unknown) bounds on x_j simply produce no tightening. Every derived
    % bound is IMPLIED by box+rows, so the result still encloses the true
    % constrained joint region (this can only shrink the box -- sound). Two-row
    % `==` pairs propagate both directions, which is what closes the g-block.
    % Monotone (bounds only shrink), so the fixed point exists; iteration is
    % capped defensively and an early stop just leaves looser (still sound)
    % bounds.
    if isempty(C), return; end
    maxIter = 8 + 2*size(C, 1);
    iter = 0;
    changed = true;
    while changed && iter < maxIter
        iter = iter + 1;
        changed = false;
        for r = 1:size(C, 1)
            nz = find(C(r, :));
            if numel(nz) ~= 2, continue; end   % only two-term rows are propagated
            for s = 1:2
                i = nz(s); j = nz(3-s);
                ai = C(r, i); aj = C(r, j);
                if aj > 0, mj = aj*lb(j); else, mj = aj*ub(j); end
                if isnan(mj), continue; end
                bnd = (d(r) - mj)/ai;
                if ai > 0
                    if isnan(ub(i)) || bnd < ub(i)
                        ub(i) = bnd; changed = true;
                    end
                else
                    if isnan(lb(i)) || bnd > lb(i)
                        lb(i) = bnd; changed = true;
                    end
                end
            end
        end
    end
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
