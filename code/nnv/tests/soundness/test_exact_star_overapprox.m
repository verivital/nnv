classdef test_exact_star_overapprox < matlab.unittest.TestCase
    % SOUNDNESS [42]: an OVER-APPROXIMATE reach can prove a property SAFE/robust
    % definitively, but can NEVER certify it UNSAFE/not-robust (the offending
    % region may exist only in the over-approximation). NNV's exact-star verdict
    % logic promotes unknown -> not-robust; that promotion is valid ONLY when the
    % whole reach was exact. These tests pin that: a network containing a layer
    % with no exact reachability (e.g. GeluLayer) must NOT return a definitive
    % not-robust under exact-star, and must emit one warning.

    methods (Test)

        function test_exact_capable_net_keeps_exact(testCase)
            % FC + ReLU is exact under exact-star -> exactReach true, and a
            % genuinely-violated property may be reported definitively (0).
            net = NN({FullyConnectedLayer('fc', eye(2), [0;0]), ReluLayer('relu')}, []);
            IS  = Star([-1;-1], [1;1]);
            ro  = struct; ro.reachMethod = 'exact-star';
            U   = HalfSpace([-1 0], -0.3);            % unsafe: y1 >= 0.3 (set straddles it)
            res = net.verify_robustness(IS, ro, U);
            testCase.verifyTrue(net.exactReach, 'FC+ReLU exact-star must be exact');
            testCase.verifyEqual(res, 0, 'exact net may report definitive not-robust');
        end

        function test_overapprox_net_cannot_be_not_robust(testCase)
            % FC + Gelu is over-approximate -> exactReach false; verify must
            % return 1 or 2 but NEVER a definitive 0 under exact-star.
            net = NN({FullyConnectedLayer('fc', eye(2), [0;0]), GeluLayer('g')}, []);
            IS  = Star([-1;-1], [1;1]);
            ro  = struct; ro.reachMethod = 'exact-star';
            U   = HalfSpace([-1 0], -0.3);
            warnState = warning('off', 'NNV:exactStarOverapprox');
            cleanup = onCleanup(@() warning(warnState));
            res = net.verify_robustness(IS, ro, U);
            testCase.verifyFalse(net.exactReach, 'FC+Gelu exact-star must be over-approximate');
            testCase.verifyNotEqual(res, 0, ...
                'SOUNDNESS: over-approximation must not certify not-robust');
            testCase.verifyEqual(res, 2, 'intersecting over-approx -> unknown');
        end

        function test_overapprox_can_still_prove_robust(testCase)
            % Over-approximation CAN prove robustness definitively (if the
            % over-approx is safe, the exact set is safe too).
            net = NN({FullyConnectedLayer('fc', eye(2), [0;0]), GeluLayer('g')}, []);
            IS  = Star([-1;-1], [1;1]);
            ro  = struct; ro.reachMethod = 'exact-star';
            Usafe = HalfSpace([1 0], -100);           % unsafe: y1 <= -100 (never reached)
            warnState = warning('off', 'NNV:exactStarOverapprox');
            cleanup = onCleanup(@() warning(warnState));
            res = net.verify_robustness(IS, ro, Usafe);
            testCase.verifyEqual(res, 1, 'over-approx must still prove robust');
        end

        function test_warns_once_on_degraded_exact_star(testCase)
            % One warning (NNV:exactStarOverapprox) is raised when exact-star is
            % requested but cannot be honoured.
            net = NN({FullyConnectedLayer('fc', eye(2), [0;0]), GeluLayer('g')}, []);
            IS  = Star([-1;-1], [1;1]);
            ro  = struct; ro.reachMethod = 'exact-star';
            testCase.verifyWarning(@() net.reach(IS, ro), 'NNV:exactStarOverapprox');
        end

        function test_approx_star_never_claims_exact(testCase)
            % approx-star is already sound (never promotes to not-robust); the
            % flag must be false and no warning raised.
            net = NN({FullyConnectedLayer('fc', eye(2), [0;0]), ReluLayer('relu')}, []);
            IS  = Star([-1;-1], [1;1]);
            ro  = struct; ro.reachMethod = 'approx-star';
            net.reach(IS, ro);
            testCase.verifyFalse(net.exactReach);
        end

        function test_layer_exactness_classification(testCase)
            % Spot-check the exact/over-approx classifier used by the gate.
            net = NN({ReluLayer('r')}, []);   % any net just to reach the method
            testCase.verifyTrue(net.layer_reach_is_exact(FullyConnectedLayer('fc', eye(2), [0;0])));
            testCase.verifyTrue(net.layer_reach_is_exact(ReluLayer('r')));
            testCase.verifyFalse(net.layer_reach_is_exact(GeluLayer('g')));
            testCase.verifyFalse(net.layer_reach_is_exact(ElementwiseProductLayer('ep')));
            testCase.verifyFalse(net.layer_reach_is_exact(SiLULayer('si')));
            % final softmax = identity (exact); intermediate softmax = over-approx
            sf = SoftmaxLayer('sf'); sf.IsFinalLayer = true;
            sm = SoftmaxLayer('sm'); sm.IsFinalLayer = false;
            testCase.verifyTrue(net.layer_reach_is_exact(sf));
            testCase.verifyFalse(net.layer_reach_is_exact(sm));
        end

    end
end
