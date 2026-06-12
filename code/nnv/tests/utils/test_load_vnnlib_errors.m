function tests = test_load_vnnlib_errors
% test_load_vnnlib_errors
% Function-based tests (so try/catch + verifyError are allowed) that PIN the
% CURRENT behavior of load_vnnlib on forms it does NOT fully support, and that
% assert well-formedness of REAL multi-region benchmark vnnlib files.
%
% These document behavior as-is; where current behavior is a limitation (not a
% deliberate design choice) it is flagged in the accompanying report, NOT fixed
% here.
%
% To run: results = runtests('test_load_vnnlib_errors')
    tests = functiontests(localfunctions);
end

% ---------------------------------------------------------------------------
% (g) top-level (assert (and ...)) with NO surrounding (or ...) currently ERRORS.
%     process_assertion hits the `startsWith(tline,'(and')` branch and throws:
%       "Property not supported for now. Not allowed -> assertion starting with
%        AND and no OR."
%     Pinned so a future change that starts SUPPORTING bare top-level AND (or that
%     changes the message) trips this test and gets a deliberate review.
% ---------------------------------------------------------------------------
function test_top_level_and_errors(testCase)
    tf = [tempname '.vnnlib'];
    fid = fopen(tf, 'w');
    fprintf(fid, '(declare-const X_0 Real)\n');
    fprintf(fid, '(declare-const Y_0 Real)\n(declare-const Y_1 Real)\n');
    fprintf(fid, '(assert (<= X_0 1.0))\n(assert (>= X_0 -1.0))\n');
    fprintf(fid, '(assert (and (>= Y_0 1.0) (<= Y_1 0.0)))\n');
    fclose(fid);
    c = onCleanup(@() delete(tf));

    threw = false;
    msg = '';
    try
        load_vnnlib(tf);
    catch ME
        threw = true;
        msg = ME.message;
    end
    verifyTrue(testCase, threw, ...
        'top-level (assert (and ...)) with no (or ...) should currently error');
    % Stable substring of the current error message.
    verifyTrue(testCase, contains(msg, 'assertion starting with AND'), ...
        sprintf('unexpected error message: %s', msg));
end

% ---------------------------------------------------------------------------
% (h) linear-combination constraint (assert (<= (+ Y_0 Y_1) 5)) is NOT supported.
%     process_constraint splits the constraint on whitespace and treats token 2
%     as a "Y_a" variable; for "(+" there is no "_a" index, so str2double indexing
%     throws. PIN as: load_vnnlib raises an error (does NOT silently mis-parse).
%     >>> FLAGGED IN REPORT: load_vnnlib cannot parse linear-combination output
%         constraints; competitions that use them would need a parser extension.
% ---------------------------------------------------------------------------
function test_linear_combination_unsupported(testCase)
    tf = [tempname '.vnnlib'];
    fid = fopen(tf, 'w');
    fprintf(fid, '(declare-const X_0 Real)\n');
    fprintf(fid, '(declare-const Y_0 Real)\n(declare-const Y_1 Real)\n');
    fprintf(fid, '(assert (<= X_0 1.0))\n(assert (>= X_0 -1.0))\n');
    fprintf(fid, '(assert (<= (+ Y_0 Y_1) 5.0))\n');
    fclose(fid);
    c = onCleanup(@() delete(tf));

    % Pin: this currently ERRORS (does not produce a well-formed HalfSpace).
    verifyError(testCase, @() load_vnnlib(tf), '', ...
        'linear-combination output constraint should currently error (not silently mis-parse)');
end

% ---------------------------------------------------------------------------
% (i1) REAL benchmark: acasxu prop_6 -- (assert (or (and X ...) (and X ...)))
%      gives MULTIPLE input regions, followed by an (or (and Y ...)) output spec.
%      Asserts the multi-input CELL shape and well-formedness.
% ---------------------------------------------------------------------------
function test_real_multiregion_acasxu_prop6(testCase)
    vf = acas_vnnlib('prop_6.vnnlib');
    verifyTrue(testCase, isfile(vf), sprintf('missing benchmark file: %s', vf));
    property = load_vnnlib(vf);

    % prop_6's input (or (and X..)(and X..)) -> lb/ub are CELL arrays (one per region).
    verifyTrue(testCase, iscell(property.lb), 'prop_6 lb should be a cell (multi-region)');
    verifyTrue(testCase, iscell(property.ub), 'prop_6 ub should be a cell (multi-region)');
    verifyEqual(testCase, numel(property.lb), 2, 'prop_6 has two input regions');
    verifyEqual(testCase, numel(property.ub), 2, 'prop_6 has two input regions');

    % Each region: 5 finite dims, lb <= ub elementwise.
    for r = 1:numel(property.lb)
        lb = property.lb{r}; ub = property.ub{r};
        verifyEqual(testCase, numel(lb), 5, 'acasxu has 5 inputs');
        verifyEqual(testCase, numel(ub), 5, 'acasxu has 5 inputs');
        verifyTrue(testCase, all(isfinite(lb)) && all(isfinite(ub)), 'bounds finite');
        verifyTrue(testCase, all(lb(:) <= ub(:)), 'lb <= ub elementwise');
    end

    % Single output spec; its Hg is an OR of disjuncts; every element is a HalfSpace.
    verifyEqual(testCase, numel(property.prop), 1, 'prop_6 has a single output spec');
    Hg = property.prop{1}.Hg;
    verifyTrue(testCase, isa(Hg, 'HalfSpace'), 'Hg is HalfSpace array');
    verifyTrue(testCase, numel(Hg) >= 1, 'at least one disjunct');
    assert_all_halfspaces(testCase, Hg);
end

% ---------------------------------------------------------------------------
% (i2) REAL benchmark: acasxu prop_7 -- output is (or (and Y..Y..Y..)(and Y..)).
%      Single input box (vector lb/ub), output OR with 2 disjuncts of 3 rows each.
% ---------------------------------------------------------------------------
function test_real_or_of_and_acasxu_prop7(testCase)
    vf = acas_vnnlib('prop_7.vnnlib');
    verifyTrue(testCase, isfile(vf), sprintf('missing benchmark file: %s', vf));
    property = load_vnnlib(vf);

    % Single input region -> plain vectors.
    verifyFalse(testCase, iscell(property.lb), 'prop_7 has a single input box');
    verifyEqual(testCase, numel(property.lb), 5, 'acasxu 5 inputs');
    verifyTrue(testCase, all(property.lb(:) <= property.ub(:)), 'lb <= ub');

    verifyEqual(testCase, numel(property.prop), 1, 'single output spec');
    Hg = property.prop{1}.Hg;
    verifyEqual(testCase, numel(Hg), 2, 'prop_7 has two OR disjuncts');
    % each disjunct is (and (<= ..)(<= ..)(<= ..)) -> 3 rows
    verifyEqual(testCase, size(Hg(1).G, 1), 3, 'disjunct 1 has 3 AND rows');
    verifyEqual(testCase, size(Hg(2).G, 1), 3, 'disjunct 2 has 3 AND rows');
    assert_all_halfspaces(testCase, Hg);
end

% ---------------------------------------------------------------------------
% (i3) REAL benchmark: acasxu prop_3 -- four standalone (<= Y_0 Y_k) asserts
%      AND-ed into a single disjunct of 4 rows. Pins the "series of asserts -> one
%      polytope" path on a real file.
% ---------------------------------------------------------------------------
function test_real_series_of_asserts_acasxu_prop3(testCase)
    vf = acas_vnnlib('prop_3.vnnlib');
    verifyTrue(testCase, isfile(vf), sprintf('missing benchmark file: %s', vf));
    property = load_vnnlib(vf);

    verifyEqual(testCase, numel(property.prop), 1, 'single output spec');
    Hg = property.prop{1}.Hg;
    verifyEqual(testCase, numel(Hg), 1, 'prop_3 output is a single conjunction (no OR)');
    verifyEqual(testCase, size(Hg(1).G, 1), 4, 'four AND-ed (<= Y_0 Y_k) rows');
    assert_all_halfspaces(testCase, Hg);
end

% ---------------------------------------------------------------------------
% Helpers
% ---------------------------------------------------------------------------
function vf = acas_vnnlib(name)
    vf = fullfile(nnvroot(), 'code', 'nnv', 'examples', 'NNV2.0', 'Submission', ...
        'CAV2023', 'NNV_vs_MATLAB', 'acas', 'vnnlib', name);
end

function assert_all_halfspaces(testCase, Hg)
    % Every element of an Hg array must be a HalfSpace with consistent G/g shapes.
    verifyTrue(testCase, isa(Hg, 'HalfSpace'), 'Hg must be a HalfSpace array');
    for i = 1:numel(Hg)
        verifyTrue(testCase, isa(Hg(i), 'HalfSpace'), 'each disjunct is a HalfSpace');
        verifyEqual(testCase, size(Hg(i).G, 1), numel(Hg(i).g), ...
            'rows of G must match length of g');
        verifyTrue(testCase, size(Hg(i).G, 1) >= 1, 'at least one constraint row');
    end
end
