function tests = test_load_vnnlib_audit_fixes
% test_load_vnnlib_audit_fixes
% Regression tests for three latent soundness bugs in load_vnnlib found by the
% VNN-COMP soundness audit (2026-06-11) and FIXED. Each test pins the CORRECT
% post-fix behavior; the pre-fix behavior is noted in the comment so a regression
% is unambiguous. All three are latent on the 2025 benchmark set (which dodges
% them) but would flip verdicts on a 2026 benchmark using the form.
%
% To run: results = runtests('test_load_vnnlib_audit_fixes')
    tests = functiontests(localfunctions);
end

% ---------------------------------------------------------------------------
% Finding 2: a STRICT '>' input constraint must set the LOWER bound.
% Pre-fix: process_input_constraint routed to lower only on '>=', so '(> X_i v)'
% fell to the else branch and set the UPPER bound -> wrong (shifted) input box ->
% the falsifier/reach ran on the wrong domain -> false sat or false unsat.
% ---------------------------------------------------------------------------
function test_strict_gt_input_is_lower_bound(testCase)
    tf = [tempname '.vnnlib']; c = onCleanup(@() delete(tf));
    fid = fopen(tf, 'w');
    fprintf(fid, '(declare-const X_0 Real)\n(declare-const Y_0 Real)\n');
    fprintf(fid, '(assert (<= X_0 1.0))\n(assert (> X_0 0.6))\n');
    fprintf(fid, '(assert (<= Y_0 0.0))\n');
    fclose(fid);
    P = load_vnnlib(tf);
    verifyEqual(testCase, double(P.lb(1)), 0.6, 'AbsTol', 1e-6, ...
        'strict (> X_0 0.6) must set the LOWER bound to 0.6');
    verifyEqual(testCase, double(P.ub(1)), 1.0, 'AbsTol', 1e-6, ...
        'the (<= X_0 1.0) upper bound must be preserved');
    verifyTrue(testCase, P.lb(1) <= P.ub(1), 'lb <= ub');
end

% ---------------------------------------------------------------------------
% Finding 5: a standalone conjunct asserted BEFORE an (or ...) must be AND-ed
% into every disjunct (mirror of the already-fixed after-the-or case).
% Pre-fix: the else branch did [last_ast.Hg.G; ast.Hg.G] with a multi-element
% ast.Hg -> "Dimensions of arrays being concatenated are not consistent" error
% (lost the instance), or a silent mis-build dropping the disjunction.
% ---------------------------------------------------------------------------
function test_conjunct_before_or_is_conjoined(testCase)
    tf = [tempname '.vnnlib']; c = onCleanup(@() delete(tf));
    fid = fopen(tf, 'w');
    fprintf(fid, '(declare-const X_0 Real)\n');
    fprintf(fid, '(declare-const Y_0 Real)\n(declare-const Y_1 Real)\n(declare-const Y_2 Real)\n');
    fprintf(fid, '(assert (<= X_0 1.0))\n(assert (>= X_0 0.0))\n');
    fprintf(fid, '(assert (<= Y_0 0.0))\n');                                 % conjunct FIRST
    fprintf(fid, '(assert (or (and (>= Y_1 1.0)) (and (>= Y_2 1.0))))\n');    % then the OR
    fclose(fid);
    P = load_vnnlib(tf);
    Hg = P.prop{1}.Hg;
    verifyEqual(testCase, numel(Hg), 2, 'two OR disjuncts');
    verifyEqual(testCase, size(Hg(1).G, 1), 2, 'disjunct 1: the Y_0<=0 conjunct AND its own row');
    verifyEqual(testCase, size(Hg(2).G, 1), 2, 'disjunct 2: the Y_0<=0 conjunct AND its own row');
    % the conjunct constrains Y_0 (index 1) -> a nonzero coefficient in column 1 of every disjunct
    verifyTrue(testCase, any(Hg(1).G(:, 1) ~= 0), 'Y_0 conjunct present in disjunct 1');
    verifyTrue(testCase, any(Hg(2).G(:, 1) ~= 0), 'Y_0 conjunct present in disjunct 2');
end

% ---------------------------------------------------------------------------
% Finding 1: in a multi-input-region (or ...), each disjunct's box must start
% from the GLOBAL declared bounds, not inherit the previous disjunct's values.
% Pre-fix: process_multiple_inputs never reset lb_input/ub_input between
% disjuncts, so a dimension constrained in region 1 LEAKED into region 2 when
% region 2 omitted it -> region 2's box was wrong -> false sat/unsat.
% ---------------------------------------------------------------------------
function test_multiinput_no_bound_leak(testCase)
    tf = [tempname '.vnnlib']; c = onCleanup(@() delete(tf));
    fid = fopen(tf, 'w');
    fprintf(fid, '(declare-const X_0 Real)\n(declare-const X_1 Real)\n(declare-const Y_0 Real)\n');
    % region 1 constrains X_0 and X_1; region 2 constrains ONLY X_0.
    fprintf(fid, ['(assert (or (and (<= X_0 1.0) (>= X_0 0.0) (<= X_1 1.0) (>= X_1 0.0)) ' ...
                  '(and (<= X_0 5.0) (>= X_0 4.0))))\n']);
    fprintf(fid, '(assert (>= Y_0 0.0))\n');
    fclose(fid);
    P = load_vnnlib(tf);
    verifyTrue(testCase, iscell(P.lb) && numel(P.lb) == 2, 'two input regions (cell lb/ub)');
    % region 1: X_1 in [0,1] as specified.
    verifyEqual(testCase, double(P.lb{1}(2)), 0.0, 'AbsTol', 1e-6, 'region 1 X_1 lb');
    verifyEqual(testCase, double(P.ub{1}(2)), 1.0, 'AbsTol', 1e-6, 'region 1 X_1 ub');
    % region 2: X_1 must NOT have leaked region 1's [0,1]; it falls back to the
    % global default (here the unset template, 0/0), so ub{2}(2) ~= 1.0.
    verifyNotEqual(testCase, double(P.ub{2}(2)), 1.0, ...
        'region 2 X_1 must NOT inherit region 1''s upper bound (no leak)');
    verifyEqual(testCase, double(P.lb{2}(1)), 4.0, 'AbsTol', 1e-6, 'region 2 X_0 lb as specified');
    verifyEqual(testCase, double(P.ub{2}(1)), 5.0, 'AbsTol', 1e-6, 'region 2 X_0 ub as specified');
end

% ---------------------------------------------------------------------------
% Finding 1b: the SAME no-leak rule for the COMBINED input/output (or ...) form
% [load_vnnlib "option 3", process_combined_input_output]. Each disjunct pairs an
% input region with an output halfspace; its input box must start from the GLOBAL
% declared bounds, not inherit the previous disjunct's. Pre-fix:
% process_combined_input_output never reset lb_input/ub_input between disjuncts, so
% an X dimension constrained in disjunct 1 LEAKED into disjunct 2 -> wrong (too
% narrow) box -> false sat/unsat (-150). Mirrors the process_multiple_inputs fix.
% ---------------------------------------------------------------------------
function test_combined_input_output_no_bound_leak(testCase)
    tf = [tempname '.vnnlib']; c = onCleanup(@() delete(tf));
    fid = fopen(tf, 'w');
    fprintf(fid, '(declare-const X_0 Real)\n(declare-const X_1 Real)\n(declare-const Y_0 Real)\n');
    fprintf(fid, '(assert (<= X_0 10.0))\n(assert (>= X_0 -10.0))\n');   % global box
    fprintf(fid, '(assert (<= X_1 10.0))\n(assert (>= X_1 -10.0))\n');
    % disjunct 1 constrains X_0 (+ a Y -> option 3); disjunct 2 constrains ONLY X_1.
    fprintf(fid, ['(assert (or (and (>= X_0 1.0) (<= X_0 2.0) (>= Y_0 0.0)) ' ...
                  '(and (>= X_1 3.0) (<= X_1 4.0) (>= Y_0 0.0))))\n']);
    fclose(fid);
    P = load_vnnlib(tf);
    verifyTrue(testCase, iscell(P.lb) && numel(P.lb) == 2, 'two input regions (option 3)');
    % disjunct 1: X_0 in [1,2] as specified
    verifyEqual(testCase, double(P.lb{1}(1)), 1.0, 'AbsTol', 1e-6, 'disjunct1 X_0 lb');
    verifyEqual(testCase, double(P.ub{1}(1)), 2.0, 'AbsTol', 1e-6, 'disjunct1 X_0 ub');
    % disjunct 2: X_0 must NOT have leaked [1,2]; it stays at the GLOBAL [-10,10].
    verifyEqual(testCase, double(P.lb{2}(1)), -10.0, 'AbsTol', 1e-6, ...
        'disjunct2 X_0 must NOT inherit disjunct1''s lower bound (no leak)');
    verifyEqual(testCase, double(P.ub{2}(1)), 10.0, 'AbsTol', 1e-6, ...
        'disjunct2 X_0 must NOT inherit disjunct1''s upper bound (no leak)');
    verifyEqual(testCase, double(P.lb{2}(2)), 3.0, 'AbsTol', 1e-6, 'disjunct2 X_1 lb as specified');
    verifyEqual(testCase, double(P.ub{2}(2)), 4.0, 'AbsTol', 1e-6, 'disjunct2 X_1 ub as specified');
end
