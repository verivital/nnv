classdef test_pgd_falsify < matlab.unittest.TestCase
    % Tests for the gradient-based falsifier (pgd_falsify) and the SAT-witness
    % validator (validate_witness) -- VNN-COMP 2026 strategy Pillars 1 & 2.
    %
    % Net under test: y = x1 + x2 (a 2->1 linear dlnetwork) on the box [0,1]^2.
    % Unsafe region: y >= 1.5, encoded as the HalfSpace {y : -y <= -1.5}. A
    % counterexample exists (e.g. x=[1;1] -> y=2), so PGD must find one, the
    % validator must confirm it, and must REJECT non-violating / out-of-box points.

    methods (Static)
        function net = linsum_net()
            layers = [ featureInputLayer(2, 'Name', 'in')
                       fullyConnectedLayer(1, 'Name', 'fc', 'Weights', [1 1], 'Bias', 0) ];
            net = dlnetwork(layers);
        end
    end

    methods (Test)

        function test_pgd_finds_counterexample(testCase)
            net = test_pgd_falsify.linsum_net();
            lb = [0;0]; ub = [1;1];
            Hs = HalfSpace(-1, -1.5);                 % {y : y >= 1.5}
            opts = struct('seed', 0, 'n_restarts', 12, 'n_steps', 30, 'lr', 0.2);
            [cex, found] = pgd_falsify(net, lb, ub, Hs, 2, 'CB', opts);
            testCase.verifyTrue(found, 'PGD should find a counterexample for y>=1.5 on [0,1]^2');
            testCase.verifyEqual(numel(cex), 2, 'cex = {x; y}');
            % the returned point is a genuine violation
            testCase.verifyTrue(validate_witness(net, cex{1}, lb, ub, Hs, 2, 'CB'), ...
                'PGD-found witness must validate');
        end

        function test_validate_rejects_non_violating(testCase)
            net = test_pgd_falsify.linsum_net();
            lb = [0;0]; ub = [1;1]; Hs = HalfSpace(-1, -1.5);
            % x=[0;0] -> y=0, NOT in {y>=1.5}
            testCase.verifyFalse(validate_witness(net, [0;0], lb, ub, Hs, 2, 'CB'), ...
                'a non-violating point must be rejected (else -150 risk)');
        end

        function test_validate_rejects_out_of_box(testCase)
            net = test_pgd_falsify.linsum_net();
            lb = [0;0]; ub = [1;1]; Hs = HalfSpace(-1, -1.5);
            % x=[2;2] violates the OUTPUT spec (y=4>=1.5) but is OUTSIDE the input box
            testCase.verifyFalse(validate_witness(net, [2;2], lb, ub, Hs, 2, 'CB'), ...
                'an out-of-box witness must be rejected');
        end

        function test_pgd_no_op_on_non_dlnetwork(testCase)
            % On a non-dlnetwork (e.g. an NNV NN), pgd must gracefully return
            % found=false so the caller falls back to random sampling.
            Hs = HalfSpace(-1, -1.5);
            [~, found] = pgd_falsify("not a dlnetwork", [0;0], [1;1], Hs, 2, 'CB', struct());
            testCase.verifyFalse(found);
        end

    end
end
