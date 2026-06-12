function tests = test_vnnlib2_multinet
%TEST_VNNLIB2_MULTINET  Phase 2 (VNN-LIB 2.0 multi-network `equal-to`) tests:
%   (a) load_vnnlib2 parses a monotonic_acasxu-style two-network spec into
%       property.multinet: joint stacked box (g-block bounds CLOSED by interval
%       propagation through the coupling rows), jointC/jointd coupling rows over
%       the CONCRETE stacked input, and the cross-network unsafe HalfSpace with
%       the single-net polarity convention (asserted comparison maps DIRECTLY
%       into G*y <= g -- no negation; see mn_linear_two_terms);
%   (b) multinet_product_net builds the weight-shared product net (blkdiag
%       weights, stacked bias, pass-through ReLU) and the product evaluates
%       EXACTLY like running f on each half;
%   (c) verify_multinet soundly refuses (status 2) non-whitelisted layers and
%       non-`equal` equivalence kinds -- the -150 guard.
%
%   Gating discipline: property.unsupported stays true unless EVERY assert
%   parsed cleanly AND the joint box closed, so under-constrained / 3-network /
%   isomorphic-to / shape-mismatched files all stay `unknown`.
%
%   Self-contained (inline temp specs, hand-built tiny nets); no reach / LP
%   solver is exercised, so this runs fast in CI.
%
%   To run: results = runtests('test_vnnlib2_multinet')
    tests = functiontests(localfunctions);
end

% ---------- (a) parser: monotonic_acasxu-style equal-to pair ----------

function test_multinet_parse_monotonic_style(tc)
    % two 5-d networks, equal-to; box on X_f; == couplings X_f[i]=X_g[i] for
    % i~=0; X_f[0] >= X_g[0] (+ lower bound on X_g[0]); output (< Y_f[3] Y_g[3])
    tf = write_vnnlib2({ ...
        '(vnnlib-version <2.0>)', ...
        '(declare-network f (declare-input X_f real [5]) (declare-output Y_f real [5]))', ...
        '(declare-network g (equal-to f) (declare-input X_g real [5]) (declare-output Y_g real [5]))', ...
        '(assert (and (<= X_f[0] 0.6) (>= X_f[0] -0.1)))', ...
        '(assert (and (<= X_f[1] 0.0) (>= X_f[1] -0.25)))', ...
        '(assert (and (<= X_f[2] 0.5) (>= X_f[2] 0.25)))', ...
        '(assert (== X_f[3] 0.2))', ...
        '(assert (== X_f[4] 0.25))', ...
        '(assert (and (>= X_f[0] X_g[0]) (>= X_g[0] -0.1)))', ...
        '(assert (== X_f[1] X_g[1]))', ...
        '(assert (== X_f[2] X_g[2]))', ...
        '(assert (== X_f[3] X_g[3]))', ...
        '(assert (== X_f[4] X_g[4]))', ...
        '(assert (< Y_f[3] Y_g[3]))'});
    property = load_vnnlib2(tf); delete(tf);
    verifyFalse(tc, is_unsupported(property), 'clean equal-to pair must parse as supported');
    verifyTrue(tc, isfield(property, 'multinet'));
    mn = property.multinet;
    verifyEqual(tc, mn.names, {'f', 'g'});
    verifyEqual(tc, mn.equivKind, 'equal');
    verifyEqual(tc, mn.inShapes{1}, 5);
    verifyEqual(tc, mn.inShapes{2}, 5);
    verifyEqual(tc, mn.outShapes{1}, 5);
    verifyEqual(tc, mn.outShapes{2}, 5);

    % joint box (10x1, f-block then g-block). The g-block has NO direct bounds
    % except lb(X_g[0]); everything else must be PROPAGATED through the
    % couplings: X_g[0] <= X_f[0] <= 0.6 and X_g[i] == X_f[i] for i = 1..4.
    expLb = [-0.1; -0.25; 0.25; 0.2; 0.25; -0.1; -0.25; 0.25; 0.2; 0.25];
    expUb = [ 0.6;  0.0;  0.5; 0.2; 0.25;  0.6;  0.0;  0.5; 0.2; 0.25];
    verifyEqual(tc, size(mn.jointLb), [10 1]);
    verifyEqual(tc, size(mn.jointUb), [10 1]);
    verifyEqual(tc, double(mn.jointLb), expLb, 'AbsTol', 1e-6);
    verifyEqual(tc, double(mn.jointUb), expUb, 'AbsTol', 1e-6);

    % coupling rows over the CONCRETE stacked x = [X_f; X_g]:
    % 1 row from (>= X_f[0] X_g[0]) + 2 rows per == coupling (4 of them) = 9
    verifyEqual(tc, size(mn.jointC), [9 10]);
    verifyEqual(tc, double(mn.jointd), zeros(9, 1), 'AbsTol', 1e-12);
    % (>= X_f[0] X_g[0])  ->  -x_f0 + x_g0 <= 0  (ONE row)
    verifyEqual(tc, double(mn.jointC(1, :)), [-1 0 0 0 0 1 0 0 0 0], 'AbsTol', 1e-12);
    % (== X_f[1] X_g[1])  ->  TWO-row inequality pair (x_f1 - x_g1 <= 0, and its negation)
    verifyEqual(tc, double(mn.jointC(2, :)), [0 1 0 0 0 0 -1 0 0 0], 'AbsTol', 1e-12);
    verifyEqual(tc, double(mn.jointC(3, :)), [0 -1 0 0 0 0 1 0 0 0], 'AbsTol', 1e-12);
    % remaining pairs follow the same two-row pattern for dims 2, 3, 4
    verifyEqual(tc, double(mn.jointC(4, :)), [0 0 1 0 0 0 0 -1 0 0], 'AbsTol', 1e-12);
    verifyEqual(tc, double(mn.jointC(8, :)), [0 0 0 0 1 0 0 0 0 -1], 'AbsTol', 1e-12);

    % cross-network UNSAFE region over the stacked [Y_f; Y_g]: VNN-LIB asserts
    % encode the sat/counterexample region ITSELF, so (< Y_f[3] Y_g[3]) maps
    % DIRECTLY (strict relaxed to <=) to  Y_f[3] - Y_g[3] <= 0:
    %   +1 at stacked col 4 (Y_f[3]), -1 at stacked col 9 (Y_g[3]), g = 0.
    % (Same convention as the single-net (<= Y[a] Y[b]) -> G(a+1)=+1,G(b+1)=-1.)
    verifyEqual(tc, numel(mn.crossProp), 1);
    expG = zeros(1, 10); expG(4) = 1; expG(9) = -1;
    verifyEqual(tc, double(mn.crossProp(1).G), expG, 'AbsTol', 1e-12, ...
        'polarity must mirror the single-net direct-mapping convention');
    verifyEqual(tc, double(mn.crossProp(1).g), 0, 'AbsTol', 1e-12);

    % legacy single-network contract fields must be unusable (runner wiring for
    % multinet is Phase 3c; nothing may mistake this for a single-net box)
    verifyTrue(tc, isempty(property.lb) && isempty(property.ub) && isempty(property.prop));
end

% ---------- (a) parser: sound gating ----------

function test_multinet_gate_isomorphic(tc)
    tf = write_vnnlib2({ ...
        '(vnnlib-version <2.0>)', ...
        '(declare-network f (declare-input X_f real [2]) (declare-output Y_f real [2]))', ...
        '(declare-network g (isomorphic-to f) (declare-input X_g real [2]) (declare-output Y_g real [2]))', ...
        '(assert (and (<= X_f[0] 1.0) (>= X_f[0] -1.0)))', ...
        '(assert (== X_f[0] X_g[0]))', ...
        '(assert (< Y_f[0] Y_g[0]))'});
    property = load_vnnlib2(tf); delete(tf);
    verifyTrue(tc, is_unsupported(property), 'isomorphic-to (different weights) must stay gated');
end

function test_multinet_gate_three_networks(tc)
    tf = write_vnnlib2({ ...
        '(vnnlib-version <2.0>)', ...
        '(declare-network f (declare-input X_f real [2]) (declare-output Y_f real [2]))', ...
        '(declare-network g (equal-to f) (declare-input X_g real [2]) (declare-output Y_g real [2]))', ...
        '(declare-network h (equal-to f) (declare-input X_h real [2]) (declare-output Y_h real [2]))', ...
        '(assert (< Y_f[0] Y_g[0]))'});
    property = load_vnnlib2(tf); delete(tf);
    verifyTrue(tc, is_unsupported(property), 'three declare-network must stay gated');
end

function test_multinet_gate_unbounded_joint_box(tc)
    % couplings + output only, NO box on X_f: propagation cannot close the
    % joint box -> unsupported stays true (the unsupported-by-default rule)
    tf = write_vnnlib2({ ...
        '(vnnlib-version <2.0>)', ...
        '(declare-network f (declare-input X_f real [5]) (declare-output Y_f real [5]))', ...
        '(declare-network g (equal-to f) (declare-input X_g real [5]) (declare-output Y_g real [5]))', ...
        '(assert (>= X_f[0] X_g[0]))', ...
        '(assert (< Y_f[3] Y_g[3]))'});
    property = load_vnnlib2(tf); delete(tf);
    verifyTrue(tc, is_unsupported(property), 'unbounded joint box must stay gated');
end

function test_multinet_gate_shape_mismatch(tc)
    % equal-to requires identical shapes (2.0 standard); [5] vs [4] must gate
    tf = write_vnnlib2({ ...
        '(vnnlib-version <2.0>)', ...
        '(declare-network f (declare-input X_f real [5]) (declare-output Y_f real [5]))', ...
        '(declare-network g (equal-to f) (declare-input X_g real [4]) (declare-output Y_g real [5]))', ...
        '(assert (< Y_f[3] Y_g[3]))'});
    property = load_vnnlib2(tf); delete(tf);
    verifyTrue(tc, is_unsupported(property), 'equal-to with mismatched shapes must stay gated');
end

% ---------- (b) product net: blkdiag construction + evaluate consistency ----------

function test_product_net_blkdiag(tc)
    [f, W1, b1, W2, b2] = tiny_fc_relu_net();
    [h, ok] = multinet_product_net(f);
    verifyTrue(tc, ok, 'FC+ReLU net must be whitelisted');
    verifyEqual(tc, numel(h.Layers), 4);
    % each FC becomes blkdiag(W, W) with stacked bias [b; b]
    verifyEqual(tc, h.Layers{1}.Weights, blkdiag(W1, W1));
    verifyEqual(tc, h.Layers{1}.Bias, [b1; b1]);
    verifyEqual(tc, h.Layers{3}.Weights, blkdiag(W2, W2));
    verifyEqual(tc, h.Layers{3}.Bias, [b2; b2]);
    % ReLU passes through unchanged (elementwise on the stacked vector)
    verifyTrue(tc, isa(h.Layers{2}, 'ReluLayer'));
    verifyTrue(tc, isa(h.Layers{4}, 'ReluLayer'));
end

function test_product_net_evaluate_consistency(tc)
    % h([x1; x2]) must equal [f(x1); f(x2)] -- the product is EXACT, not an
    % approximation (fixed deterministic probe points, incl. ReLU sign flips)
    f = tiny_fc_relu_net();
    h = multinet_product_net(f);
    probes = { [ 1.0; -1.0], [ 0.3;  2.0]; ...
               [-0.5; -0.25], [ 2.0;  0.1]; ...
               [ 0.0;  0.0], [-1.5;  1.5]; ...
               [ 0.7;  0.7], [ 0.7;  0.7] };
    for p = 1:size(probes, 1)
        x1 = probes{p, 1}; x2 = probes{p, 2};
        y = h.evaluate([x1; x2]);
        yref = [f.evaluate(x1); f.evaluate(x2)];
        verifyEqual(tc, y, yref, 'AbsTol', 1e-12, ...
            'product net must replicate f on each half of the stacked input');
    end
end

function test_product_net_rejects_non_whitelisted(tc)
    net = NN({FullyConnectedLayer('fc1', eye(2), zeros(2, 1)), TanhLayer()});
    [~, ok, reason] = multinet_product_net(net);
    verifyFalse(tc, ok, 'TanhLayer must not be whitelisted');
    verifyTrue(tc, contains(reason, 'TanhLayer'));
end

% ---------- (c) verify_multinet: sound refusals ----------

function test_verify_multinet_whitelist_gate(tc)
    % a net with a non-whitelisted layer must return unknown (2) IMMEDIATELY
    % (whitelist runs before falsification and reach)
    net = NN({FullyConnectedLayer('fc1', eye(2), zeros(2, 1)), TanhLayer()});
    property = hand_built_multinet_property('equal');
    [status, counterEx] = verify_multinet(net, property, struct('reachMethod', 'approx-star'));
    verifyEqual(tc, status, 2, 'non-whitelisted layer (TanhLayer) must yield unknown');
    verifyTrue(tc, isempty(counterEx));
end

function test_verify_multinet_equivkind_gate(tc)
    % equivKind ~= 'equal' (e.g. isomorphic-to) must NEVER reach the equal-to
    % falsifier/product path (g would be a DIFFERENT network) -> unknown
    net = NN({FullyConnectedLayer('fc1', eye(2), zeros(2, 1)), ReluLayer('relu1')});
    property = hand_built_multinet_property('isomorphic');
    [status, counterEx] = verify_multinet(net, property, struct('reachMethod', 'approx-star'));
    verifyEqual(tc, status, 2, 'equivKind isomorphic must yield unknown');
    verifyTrue(tc, isempty(counterEx));
end

function test_verify_multinet_missing_multinet_gate(tc)
    % a property without the .multinet payload (e.g. single-net or gated parse)
    % must yield unknown, never a verdict
    net = NN({FullyConnectedLayer('fc1', eye(2), zeros(2, 1)), ReluLayer('relu1')});
    property = struct('unsupported', false, 'lb', [], 'ub', [], 'prop', {{}});
    [status, counterEx] = verify_multinet(net, property, struct('reachMethod', 'approx-star'));
    verifyEqual(tc, status, 2);
    verifyTrue(tc, isempty(counterEx));
end

% ---------- helpers ----------

function [f, W1, b1, W2, b2] = tiny_fc_relu_net()
    % tiny 2-2-2 FC+ReLU network built by hand
    W1 = [1 -2; 0.5 1]; b1 = [0.1; -0.2];
    W2 = [1 1; -1 2];   b2 = [0; 0.3];
    f = NN({FullyConnectedLayer('fc1', W1, b1), ReluLayer('relu1'), ...
            FullyConnectedLayer('fc2', W2, b2), ReluLayer('relu2')});
end

function property = hand_built_multinet_property(equivKind)
    % minimal well-formed property.multinet for two 2-d nets (joint dim 4)
    mn = struct();
    mn.names = {'f', 'g'};
    mn.inShapes = {2, 2};
    mn.outShapes = {2, 2};
    mn.equivKind = equivKind;
    mn.jointLb = -ones(4, 1);
    mn.jointUb = ones(4, 1);
    mn.jointC = zeros(0, 4);
    mn.jointd = zeros(0, 1);
    mn.crossProp = HalfSpace([1 0 -1 0], 0);
    property = struct();
    property.lb = []; property.ub = []; property.prop = {};
    property.unsupported = false;
    property.reason = '';
    property.multinet = mn;
end

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
