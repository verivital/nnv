classdef test_pgd_falsify < matlab.unittest.TestCase
    % Tests for the gradient-based falsifier (pgd_falsify) and the SAT-witness
    % validator (validate_witness) -- VNN-COMP 2026 strategy Pillars 1 & 2 -- across
    % feature, non-permuted image, and PERMUTED image (needReshape) inputs.

    properties (Access = private)
        AddedPath = ''   % runner dir added by TestClassSetup (removed in teardown)
    end

    methods (TestClassSetup)
        function addRunnerPath(tc)
            % pgd_falsify/validate_witness live in examples/Submission/VNN_COMP2026.
            % The matrix CI splits tests into shards and nothing guarantees that dir
            % is on the path (same latent flake test_validate_witness_onnx had), so
            % add it here to be self-sufficient.
            here = fileparts(mfilename('fullpath'));
            sub = fullfile(here, '..', '..', 'examples', 'Submission', 'VNN_COMP2026');
            if isfolder(sub) && ~contains(path, sub)
                addpath(sub);
                tc.AddedPath = sub;
            end
        end
    end

    methods (TestClassTeardown)
        function removeRunnerPath(tc)
            if ~isempty(tc.AddedPath)
                rmpath(tc.AddedPath);
            end
        end
    end

    methods (Static)
        function net = linsum_net()                  % y = x1 + x2 (feature, 2->1)
            layers = [ featureInputLayer(2, 'Name', 'in')
                       fullyConnectedLayer(1, 'Name', 'fc', 'Weights', [1 1], 'Bias', 0) ];
            net = dlnetwork(layers);
        end
        function net = img_sum_net()                 % 1x1x2 image -> y = c1 + c2
            layers = [ imageInputLayer([1 1 2], 'Name', 'in', 'Normalization', 'none')
                       fullyConnectedLayer(1, 'Name', 'fc', 'Weights', [1 1], 'Bias', 0) ];
            net = dlnetwork(layers);
        end
        function net = img_pos_net()                 % 1x2x1 image -> y = p1 + 10*p2
            % input layer is the MATLAB (post-permute) size [1 2 1]; the ONNX
            % (pre-permute) size is [2 1 1] with needReshape=1.
            layers = [ imageInputLayer([1 2 1], 'Name', 'in', 'Normalization', 'none')
                       fullyConnectedLayer(1, 'Name', 'fc', 'Weights', [1 10], 'Bias', 0) ];
            net = dlnetwork(layers);
        end
        function net = sum64_net()                   % y = sum(x), 64 inputs (>60 -> SPSA)
            layers = [ featureInputLayer(64, 'Name', 'in')
                       fullyConnectedLayer(1, 'Name', 'fc', 'Weights', ones(1,64), 'Bias', 0) ];
            net = dlnetwork(layers);
        end
    end

    methods (Test)

        function test_pgd_feature(testCase)
            net = test_pgd_falsify.linsum_net();
            lb = [0;0]; ub = [1;1]; Hs = HalfSpace(-1, -1.5);   % unsafe y >= 1.5
            opts = struct('seed', 0, 'n_restarts', 12, 'n_steps', 30, 'lr', 0.2);
            [cex, found] = pgd_falsify(net, lb, ub, Hs, 2, 'CB', 0, opts);
            testCase.verifyTrue(found, 'PGD should find a CE for y>=1.5 on [0,1]^2');
            testCase.verifyTrue(validate_witness(net, cex{1}, lb, ub, Hs, 2, 'CB', 0), ...
                'PGD-found witness must validate');
        end

        function test_pgd_image_no_permute(testCase)
            net = test_pgd_falsify.img_sum_net();
            lb = [0;0]; ub = [1;1]; Hs = HalfSpace(-1, -1.5);
            opts = struct('seed', 0, 'n_restarts', 12, 'n_steps', 30, 'lr', 0.2);
            [cex, found] = pgd_falsify(net, lb, ub, Hs, [1 1 2], 'SSCB', 0, opts);
            testCase.verifyTrue(found, 'PGD should find a CE on the 1x1x2 image net');
            testCase.verifyTrue(validate_witness(net, cex{1}, lb, ub, Hs, [1 1 2], 'SSCB', 0), ...
                'image witness must validate');
        end

        function test_pgd_image_needreshape1(testCase)
            % Permuted image (needReshape=1): the flat ONNX input must be reshaped to
            % the ONNX size [2 1 1] then permuted [2 1 3] to the MATLAB net input
            % [1 2 1]. Position-sensitive weights ([1 10]) make the permute MATTER, so
            % a wrong mapping would either error or fail the INDEPENDENT oracle below.
            net = test_pgd_falsify.img_pos_net();
            lb = [0;0]; ub = [1;1]; Hs = HalfSpace(-1, -6);     % unsafe y >= 6
            opts = struct('seed', 0, 'n_restarts', 16, 'n_steps', 40, 'lr', 0.25);
            [cex, found] = pgd_falsify(net, lb, ub, Hs, [2 1 1], 'default', 1, opts);
            testCase.verifyTrue(found, 'PGD should find a CE on the needReshape=1 image net');
            w = cex{1};
            % INDEPENDENT oracle: with the CORRECT permute the net computes y = w1 + 10*w2.
            testCase.verifyGreaterThanOrEqual(w(1) + 10*w(2), 6 - 1e-3, ...
                'witness must satisfy the true (correct-permute) spec y = w1 + 10*w2 >= 6');
            testCase.verifyTrue(validate_witness(net, w, lb, ub, Hs, [2 1 1], 'default', 1), ...
                'needReshape=1 witness must validate (round-trip consistent)');
        end

        function test_validate_rejects_non_violating(testCase)
            net = test_pgd_falsify.linsum_net();
            lb = [0;0]; ub = [1;1]; Hs = HalfSpace(-1, -1.5);
            testCase.verifyFalse(validate_witness(net, [0;0], lb, ub, Hs, 2, 'CB', 0), ...
                'a non-violating point must be rejected (else -150 risk)');
        end

        function test_validate_rejects_out_of_box(testCase)
            net = test_pgd_falsify.linsum_net();
            lb = [0;0]; ub = [1;1]; Hs = HalfSpace(-1, -1.5);
            testCase.verifyFalse(validate_witness(net, [2;2], lb, ub, Hs, 2, 'CB', 0), ...
                'an out-of-box witness must be rejected');
        end

        function test_pgd_no_op_on_non_dlnetwork(testCase)
            Hs = HalfSpace(-1, -1.5);
            [~, found] = pgd_falsify("not a dlnetwork", [0;0], [1;1], Hs, 2, 'CB', 0, struct());
            testCase.verifyFalse(found);
        end

        % ---- NN-manifest path: numerical-gradient PGD (no autodiff) ----

        function test_pgd_nn_finds_ce(testCase)
            % matlab2nnv -> NNV NN (the manifest path's net type); pgd_falsify must use
            % the finite-difference numerical gradient (nIn=2 <= 60) to find the CE.
            net = matlab2nnv(test_pgd_falsify.linsum_net());   % y = x1 + x2 as an NNV NN
            lb = [0;0]; ub = [1;1]; Hs = HalfSpace(-1, -1.5);  % unsafe y >= 1.5
            opts = struct('seed', 0, 'n_restarts', 12, 'n_steps', 30, 'lr', 0.2);
            [cex, found] = pgd_falsify(net, lb, ub, Hs, 2, 'default', 0, opts);
            testCase.verifyTrue(found, 'numerical PGD should find a CE for y>=1.5 on an NN net');
            testCase.verifyTrue(validate_witness(net, cex{1}, lb, ub, Hs, 2, 'default', 0), ...
                'NN-found witness must validate');
        end

        function test_pgd_nn_no_ce(testCase)
            net = matlab2nnv(test_pgd_falsify.linsum_net());
            lb = [0;0]; ub = [1;1]; Hs = HalfSpace(-1, -100);  % unsafe y >= 100: impossible
            opts = struct('seed', 0, 'n_restarts', 8, 'n_steps', 30);
            [~, found] = pgd_falsify(net, lb, ub, Hs, 2, 'default', 0, opts);
            testCase.verifyFalse(found, 'no CE exists for y>=100 on [0,1]^2');
        end

        function test_pgd_nn_spsa_highdim(testCase)
            % >60 inputs triggers the SPSA estimator (2 evals/step). y = sum(x); unsafe
            % y >= 50 on [0,1]^64 is satisfiable (x=1 -> y=64), so SPSA must climb to it.
            net = matlab2nnv(test_pgd_falsify.sum64_net());
            lb = zeros(64,1); ub = ones(64,1); Hs = HalfSpace(-1, -50);
            opts = struct('seed', 0, 'n_restarts', 6, 'n_steps', 60, 'lr', 0.3);
            [cex, found] = pgd_falsify(net, lb, ub, Hs, 64, 'default', 0, opts);
            testCase.verifyTrue(found, 'SPSA numerical PGD should find a CE for sum(x)>=50 on [0,1]^64');
            testCase.verifyTrue(validate_witness(net, cex{1}, lb, ub, Hs, 64, 'default', 0));
        end

    end
end
