function tests = test_load_vnnlib2
%TEST_LOAD_VNNLIB2  Regression tests for the VNN-LIB 2.0 parser (load_vnnlib2.m)
%   and the version-sniff dispatch added to load_vnnlib.m.
%
%   Pins BOTH the supported single-network forms (box input + linear/disjunctive
%   output, multi-dim ROW-MAJOR index flattening, ==) AND the SOUND GATING of the
%   constructs NNV cannot verify (multi-network, nonlinear, multimodal), which must
%   return property.unsupported=true so the runner emits `unknown` -- never a -150
%   wrong verdict.
%
%   Output-structure contract is identical to load_vnnlib (see test_load_vnnlib_forms):
%     property.prop{n}.Hg = array of HalfSpace; G*y <= g is the UNSAFE region; an
%     ARRAY (numel>1) is an OR; rows within one HalfSpace are an AND.
%   Sign convention (must match process_constraint):
%     (<= Y[a] c)    -> G(a+1)=+1,           g=c
%     (>= Y[a] c)    -> G(a+1)=-1,           g=-c
%     (<= Y[a] Y[b]) -> G(a+1)=+1,G(b+1)=-1, g=0
%
%   Self-contained: every spec is written to a temp file inline (no benchmark clone
%   needed), so this runs in CI. load_vnnlib2 lives in engine/utils (always on the
%   NNV path), so no setupOnce is required.
%
%   To run: results = runtests('test_load_vnnlib2')
    tests = functiontests(localfunctions);
end

% ---------- supported single-network forms ----------

function test_single_net_linear(tc)
    % box input + (<= Y[0] 0) -> 1 disjunct, 1 row
    tf = write_vnnlib2({ ...
        '(vnnlib-version <2.0>)', ...
        '(declare-network N (declare-input X float32 [2]) (declare-output Y float32 [3]))', ...
        '(assert (<= X[0] 1.0))', '(assert (>= X[0] -1.0))', ...
        '(assert (<= X[1] 2.0))', '(assert (>= X[1] -2.0))', ...
        '(assert (<= Y[0] 0.0))'});
    property = load_vnnlib2(tf); delete(tf);
    verifyFalse(tc, is_unsupported(property), 'single-net linear should be supported');
    verifyEqual(tc, double(property.lb(:)'), [-1 -2], 'AbsTol', 1e-6);
    verifyEqual(tc, double(property.ub(:)'), [ 1  2], 'AbsTol', 1e-6);
    verifyEqual(tc, numel(property.prop), 1);
    Hg = property.prop{1}.Hg;
    verifyEqual(tc, numel(Hg), 1, 'one HalfSpace, no OR');
    verifyEqual(tc, double(Hg(1).G), [1 0 0], 'AbsTol', 1e-6);
    verifyEqual(tc, double(Hg(1).g), 0, 'AbsTol', 1e-6);
end

function test_multidim_rowmajor_flatten(tc)
    % X[0,0,0,3] in shape [1,1,1,5] -> flat 3 (col 4); Y[0,2] in [1,5] -> flat 2 (col 3).
    tf = write_vnnlib2({ ...
        '(vnnlib-version <2.0>)', ...
        '(declare-network N (declare-input X float32 [1, 1, 1, 5]) (declare-output Y float32 [1, 5]))', ...
        '(assert (<= X[0,0,0,0] 0.1))', '(assert (>= X[0,0,0,0] 0.0))', ...
        '(assert (<= X[0,0,0,1] 0.2))', '(assert (>= X[0,0,0,1] 0.1))', ...
        '(assert (<= X[0,0,0,2] 0.3))', '(assert (>= X[0,0,0,2] 0.2))', ...
        '(assert (<= X[0,0,0,3] 0.5))', '(assert (>= X[0,0,0,3] 0.45))', ...
        '(assert (<= X[0,0,0,4] 0.6))', '(assert (>= X[0,0,0,4] 0.5))', ...
        '(assert (>= Y[0,2] 1.0))'});
    property = load_vnnlib2(tf); delete(tf);
    verifyFalse(tc, is_unsupported(property));
    verifyEqual(tc, numel(property.lb), 5, '5 input dims = prod(shape)');
    verifyEqual(tc, double(property.lb(4)), 0.45, 'AbsTol', 1e-6);   % flat 3 -> idx 4
    verifyEqual(tc, double(property.ub(4)), 0.50, 'AbsTol', 1e-6);
    Hg = property.prop{1}.Hg;
    verifyEqual(tc, size(Hg(1).G, 2), 5, '5 output dims');
    verifyEqual(tc, double(Hg(1).G(1,3)), -1, 'AbsTol', 1e-6);       % Y[0,2] -> col 3
    verifyEqual(tc, double(Hg(1).g(1)), -1, 'AbsTol', 1e-6);
    verifyEqual(tc, nnz(Hg(1).G(1,:)), 1, 'only Y[0,2] nonzero');
end

function test_output_or(tc)
    tf = write_vnnlib2({ ...
        '(vnnlib-version <2.0>)', ...
        '(declare-network N (declare-input X float32 [1]) (declare-output Y float32 [2]))', ...
        '(assert (<= X[0] 1.0))', '(assert (>= X[0] -1.0))', ...
        '(assert (or (and (>= Y[0] 1.0)) (and (>= Y[1] 2.0))))'});
    property = load_vnnlib2(tf); delete(tf);
    Hg = property.prop{1}.Hg;
    verifyEqual(tc, numel(Hg), 2, 'two OR disjuncts');
    verifyEqual(tc, double(Hg(1).G), [-1 0], 'AbsTol', 1e-6);
    verifyEqual(tc, double(Hg(1).g), -1, 'AbsTol', 1e-6);
    verifyEqual(tc, double(Hg(2).G), [0 -1], 'AbsTol', 1e-6);
    verifyEqual(tc, double(Hg(2).g), -2, 'AbsTol', 1e-6);
end

function test_input_equality_degenerate_box(tc)
    tf = write_vnnlib2({ ...
        '(vnnlib-version <2.0>)', ...
        '(declare-network N (declare-input X float32 [2]) (declare-output Y float32 [1]))', ...
        '(assert (== X[0] 0.5))', ...
        '(assert (<= X[1] 1.0))', '(assert (>= X[1] -1.0))', ...
        '(assert (<= Y[0] 0.0))'});
    property = load_vnnlib2(tf); delete(tf);
    verifyEqual(tc, double(property.lb(1)), 0.5, 'AbsTol', 1e-6);
    verifyEqual(tc, double(property.ub(1)), 0.5, 'AbsTol', 1e-6);
end

function test_output_equality_two_rows(tc)
    tf = write_vnnlib2({ ...
        '(vnnlib-version <2.0>)', ...
        '(declare-network N (declare-input X float32 [1]) (declare-output Y float32 [1]))', ...
        '(assert (<= X[0] 1.0))', '(assert (>= X[0] -1.0))', ...
        '(assert (== Y[0] 2.0))'});
    property = load_vnnlib2(tf); delete(tf);
    Hg = property.prop{1}.Hg;
    verifyEqual(tc, size(Hg(1).G, 1), 2, '== output -> 2 rows');
    % the two rows must encode {Y<=2} (G=+1,g=2) AND {Y>=2} (G=-1,g=-2), in either
    % order -> their unsafe region is exactly {Y=2}. Assert order-independently.
    rows = sortrows([double(Hg(1).G) double(Hg(1).g)]);
    verifyEqual(tc, rows, [-1 -2; 1 2], 'AbsTol', 1e-6, '== -> {Y<=2, Y>=2}');
end

function test_var_var_output(tc)
    tf = write_vnnlib2({ ...
        '(vnnlib-version <2.0>)', ...
        '(declare-network N (declare-input X float32 [1]) (declare-output Y float32 [2]))', ...
        '(assert (<= X[0] 1.0))', '(assert (>= X[0] -1.0))', ...
        '(assert (<= Y[0] Y[1]))'});
    property = load_vnnlib2(tf); delete(tf);
    Hg = property.prop{1}.Hg;
    verifyEqual(tc, double(Hg(1).G(1,1)), 1, 'AbsTol', 1e-6);
    verifyEqual(tc, double(Hg(1).G(1,2)), -1, 'AbsTol', 1e-6);
    verifyEqual(tc, double(Hg(1).g(1)), 0, 'AbsTol', 1e-6);
end

% ---------- sound gating (the -150 guard) ----------

function test_gate_multi_network(tc)
    tf = write_vnnlib2({ ...
        '(vnnlib-version <2.0>)', ...
        '(declare-network f (declare-input X_f real [5]) (declare-output Y_f real [5]))', ...
        '(declare-network g (equal-to f) (declare-input X_g real [5]) (declare-output Y_g real [5]))', ...
        '(assert (>= X_f[0] X_g[0]))', '(assert (< Y_f[3] Y_g[3]))'});
    property = load_vnnlib2(tf); delete(tf);
    verifyTrue(tc, is_unsupported(property), 'multi-network must be gated');
end

function test_gate_nonlinear(tc)
    tf = write_vnnlib2({ ...
        '(vnnlib-version <2.0>)', ...
        '(declare-network N (declare-input X float32 [2]) (declare-output Y float32 [1]))', ...
        '(assert (>= (* X[0] X[0]) 1.0))', ...
        '(assert (<= Y[0] 0.0))'});
    property = load_vnnlib2(tf); delete(tf);
    verifyTrue(tc, is_unsupported(property), 'nonlinear (* X X) must be gated, never linearized');
end

function test_gate_multimodal(tc)
    tf = write_vnnlib2({ ...
        '(vnnlib-version <2.0>)', ...
        '(declare-network N (declare-input X1 float32 [2]) (declare-input X2 float32 [2]) (declare-output Y float32 [1]))', ...
        '(assert (<= X1[0] 1.0))', '(assert (>= X1[0] 0.0))', ...
        '(assert (> Y[0] 0.5))'});
    property = load_vnnlib2(tf); delete(tf);
    verifyTrue(tc, is_unsupported(property), 'multimodal (2 input tensors) must be gated');
end

% ---------- dispatch ----------

function test_dispatch_from_load_vnnlib(tc)
    % load_vnnlib() must sniff the 2.0 header and route to load_vnnlib2, producing
    % an identical result to calling load_vnnlib2 directly.
    tf = write_vnnlib2({ ...
        '(vnnlib-version <2.0>)', ...
        '(declare-network N (declare-input X float32 [2]) (declare-output Y float32 [3]))', ...
        '(assert (<= X[0] 1.0))', '(assert (>= X[0] -1.0))', ...
        '(assert (<= X[1] 2.0))', '(assert (>= X[1] -2.0))', ...
        '(assert (<= Y[0] 0.0))'});
    viaDispatch = load_vnnlib(tf);
    viaDirect   = load_vnnlib2(tf);
    delete(tf);
    verifyFalse(tc, is_unsupported(viaDispatch), 'dispatch parsed the 2.0 file');
    verifyEqual(tc, double(viaDispatch.lb), double(viaDirect.lb), 'AbsTol', 1e-6);
    verifyEqual(tc, double(viaDispatch.prop{1}.Hg(1).G), double(viaDirect.prop{1}.Hg(1).G), 'AbsTol', 1e-6);
end

% ---------- output OR/conjunct distribution (must stay a SINGLE prop cell) ----------

function test_conjunct_then_or_single_cell(tc)
    % A standalone conjunct asserted BEFORE an (or ...) must be AND-ed into EVERY
    % disjunct, leaving ONE prop cell. (A second cell would be silently ignored by the
    % runner's length(prop)==1 verify path -> -150.)
    tf = write_vnnlib2({ ...
        '(vnnlib-version <2.0>)', ...
        '(declare-network N (declare-input X float32 [1]) (declare-output Y float32 [2]))', ...
        '(assert (<= X[0] 1.0))', '(assert (>= X[0] -1.0))', ...
        '(assert (>= Y[0] 0.0))', ...
        '(assert (or (and (<= Y[0] 1.0)) (and (>= Y[1] 2.0))))'});
    property = load_vnnlib2(tf); delete(tf);
    verifyFalse(tc, is_unsupported(property));
    verifyEqual(tc, numel(property.prop), 1, 'must be ONE prop cell (conjunct distributed into OR)');
    Hg = property.prop{1}.Hg;
    verifyEqual(tc, numel(Hg), 2, 'two OR disjuncts');
    verifyTrue(tc, has_row(Hg(1), [-1 0], 0), 'disjunct1 carries the conjunct Y0>=0');
    verifyTrue(tc, has_row(Hg(2), [-1 0], 0), 'disjunct2 carries the conjunct Y0>=0');
    verifyTrue(tc, has_row(Hg(1), [1 0], 1),  'disjunct1 keeps Y0<=1');
    verifyTrue(tc, has_row(Hg(2), [0 -1], -2),'disjunct2 keeps Y1>=2');
end

function test_or_then_conjunct_single_cell(tc)
    % Mirror: standalone conjunct AFTER the (or ...) -> still one cell, in every disjunct.
    tf = write_vnnlib2({ ...
        '(vnnlib-version <2.0>)', ...
        '(declare-network N (declare-input X float32 [1]) (declare-output Y float32 [2]))', ...
        '(assert (<= X[0] 1.0))', '(assert (>= X[0] -1.0))', ...
        '(assert (or (and (<= Y[0] 1.0)) (and (>= Y[1] 2.0))))', ...
        '(assert (>= Y[0] 0.0))'});
    property = load_vnnlib2(tf); delete(tf);
    verifyEqual(tc, numel(property.prop), 1, 'must be ONE prop cell');
    Hg = property.prop{1}.Hg;
    verifyEqual(tc, numel(Hg), 2, 'two OR disjuncts');
    verifyTrue(tc, has_row(Hg(1), [-1 0], 0), 'disjunct1 carries the trailing conjunct');
    verifyTrue(tc, has_row(Hg(2), [-1 0], 0), 'disjunct2 carries the trailing conjunct');
end

function test_gate_product_of_sums(tc)
    % Two separate (or ...) assertions = (OR) AND (OR) = product-of-sums, which the
    % single-disjunctive-spec contract cannot represent -> must gate to unsupported.
    tf = write_vnnlib2({ ...
        '(vnnlib-version <2.0>)', ...
        '(declare-network N (declare-input X float32 [1]) (declare-output Y float32 [2]))', ...
        '(assert (<= X[0] 1.0))', '(assert (>= X[0] -1.0))', ...
        '(assert (or (and (<= Y[0] 1.0)) (and (>= Y[0] 2.0))))', ...
        '(assert (or (and (<= Y[1] 1.0)) (and (>= Y[1] 2.0))))'});
    property = load_vnnlib2(tf); delete(tf);
    verifyTrue(tc, is_unsupported(property), 'product-of-sums must be gated');
end

% ---------- helpers ----------

function tf = write_vnnlib2(lines)
    tf = [tempname '.vnnlib'];
    fid = fopen(tf, 'w');
    for k = 1:numel(lines)
        fprintf(fid, '%s\n', lines{k});
    end
    fclose(fid);
end

function u = is_unsupported(p)
    u = isfield(p, 'unsupported') && p.unsupported;
end

function tf = has_row(hs, grow, gval)
    % true if HalfSpace hs has a row equal to grow with rhs gval (order-independent)
    G = double(hs.G); g = double(hs.g);
    tf = false;
    for r = 1:size(G, 1)
        if isequal(G(r, :), grow) && abs(g(r) - gval) < 1e-6
            tf = true; return;
        end
    end
end
