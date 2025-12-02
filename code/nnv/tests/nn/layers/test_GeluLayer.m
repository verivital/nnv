%% Test Suite for GeluLayer
% Tests for the GELU (Gaussian Error Linear Unit) activation layer
%
% Author: NNV Team
% Date: December 2025

classdef test_GeluLayer < matlab.unittest.TestCase

    methods (Test)

        %% Constructor Tests

        function test_constructor_default(testCase)
            % Test default constructor
            layer = GeluLayer();

            testCase.verifyClass(layer, 'GeluLayer');
            testCase.verifyEqual(layer.Name, 'gelu_layer');
            testCase.verifyEqual(layer.NumInputs, 1);
            testCase.verifyEqual(layer.NumOutputs, 1);
        end

        function test_constructor_with_name(testCase)
            % Test constructor with name
            layer = GeluLayer('my_gelu');

            testCase.verifyEqual(layer.Name, 'my_gelu');
        end

        function test_constructor_full(testCase)
            % Test constructor with all arguments
            layer = GeluLayer('gelu', 1, {'in'}, 1, {'out'});

            testCase.verifyEqual(layer.Name, 'gelu');
            testCase.verifyEqual(layer.NumInputs, 1);
            testCase.verifyEqual(layer.InputNames, {'in'});
        end

        %% Evaluation Tests

        function test_evaluate_zeros(testCase)
            % GELU(0) = 0
            layer = GeluLayer();
            y = layer.evaluate(0);
            testCase.verifyEqual(y, 0, 'AbsTol', 1e-10);
        end

        function test_evaluate_positive(testCase)
            % GELU(x) ≈ x for large positive x
            layer = GeluLayer();
            y = layer.evaluate(3);
            testCase.verifyGreaterThan(y, 2.9);
            testCase.verifyLessThan(y, 3.1);
        end

        function test_evaluate_negative(testCase)
            % GELU(x) → 0 for large negative x
            layer = GeluLayer();
            y = layer.evaluate(-3);
            testCase.verifyLessThan(abs(y), 0.05);
        end

        function test_evaluate_small_negative(testCase)
            % GELU has a minimum around x ≈ -0.75
            layer = GeluLayer();
            y_min = layer.evaluate(-0.75);
            testCase.verifyLessThan(y_min, 0);  % Should be negative
            testCase.verifyGreaterThan(y_min, -0.2);  % But not too negative
        end

        function test_evaluate_vector(testCase)
            % Test vector input
            layer = GeluLayer();
            x = [-2, -1, 0, 1, 2];
            y = layer.evaluate(x);

            testCase.verifySize(y, [1, 5]);
            testCase.verifyEqual(y(3), 0, 'AbsTol', 1e-10);  % GELU(0) = 0
            testCase.verifyGreaterThan(y(4), 0);  % GELU(1) > 0
            testCase.verifyGreaterThan(y(5), y(4));  % GELU(2) > GELU(1)
        end

        function test_evaluate_matrix(testCase)
            % Test matrix input
            layer = GeluLayer();
            x = randn(3, 4);
            y = layer.evaluate(x);

            testCase.verifySize(y, [3, 4]);
            testCase.verifyFalse(any(isnan(y(:))));
        end

        function test_evaluate_3d_array(testCase)
            % Test 3D array input (like CNN feature maps)
            layer = GeluLayer();
            x = randn(4, 4, 3);
            y = layer.evaluate(x);

            testCase.verifySize(y, [4, 4, 3]);
        end

        function test_evaluate_monotonicity(testCase)
            % Test monotonicity for x > -0.5
            layer = GeluLayer();
            x = linspace(0, 5, 100);
            y = layer.evaluate(x);

            dy = diff(y);
            testCase.verifyTrue(all(dy > 0), 'GELU should be increasing for x > 0');
        end

        %% Reachability Tests - Star

        function test_reach_star_positive(testCase)
            % Test reachability with positive bounds
            layer = GeluLayer();

            lb = ones(4, 1);
            ub = 2 * ones(4, 1);
            S = Star(lb, ub);

            S_out = layer.reach(S, 'approx-star');

            testCase.verifyClass(S_out, 'Star');
            testCase.verifyEqual(S_out.dim, 4);
        end

        function test_reach_star_mixed(testCase)
            % Test reachability with mixed positive/negative bounds
            layer = GeluLayer();

            lb = -ones(4, 1);
            ub = ones(4, 1);
            S = Star(lb, ub);

            S_out = layer.reach(S, 'approx-star');

            testCase.verifyClass(S_out, 'Star');
        end

        function test_reach_star_negative(testCase)
            % Test reachability with negative bounds
            layer = GeluLayer();

            lb = -2 * ones(4, 1);
            ub = -1 * ones(4, 1);
            S = Star(lb, ub);

            S_out = layer.reach(S, 'approx-star');

            testCase.verifyClass(S_out, 'Star');
        end

        %% Reachability Tests - ImageStar

        function test_reach_imagestar_bounds(testCase)
            % Test reachability with ImageStar bounds representation
            layer = GeluLayer();

            lb = zeros(4, 4);
            ub = ones(4, 4);
            IS = ImageStar(lb, ub);

            IS_out = layer.reach(IS, 'approx-star');

            testCase.verifyClass(IS_out, 'ImageStar');
            testCase.verifyFalse(isempty(IS_out.im_lb));
            testCase.verifyFalse(isempty(IS_out.im_ub));
        end

        function test_reach_imagestar_3d(testCase)
            % Test reachability with 3D ImageStar
            layer = GeluLayer();

            lb = -ones(4, 4, 3);
            ub = ones(4, 4, 3);
            IS = ImageStar(lb, ub);

            IS_out = layer.reach(IS, 'approx-star');

            testCase.verifyClass(IS_out, 'ImageStar');
            testCase.verifySize(IS_out.im_lb, [4, 4, 3]);
        end

        %% Reachability Tests - Zonotope

        function test_reach_zono(testCase)
            % Test zonotope reachability
            layer = GeluLayer();

            center = zeros(4, 1);
            generators = eye(4);
            Z = Zono(center, generators);

            Z_out = layer.reach(Z, 'approx-zono');

            testCase.verifyClass(Z_out, 'Zono');
        end

        function test_reach_imagezono(testCase)
            % Test ImageZono reachability
            layer = GeluLayer();

            lb = -ones(4, 4);
            ub = ones(4, 4);
            IZ = ImageZono(lb, ub);

            IZ_out = layer.reach(IZ, 'approx-zono');

            testCase.verifyClass(IZ_out, 'ImageZono');
        end

        %% NN Calling Convention Tests

        function test_reach_nn_convention_7args(testCase)
            % Test reach with 7 arguments (NN.reach calling convention)
            layer = GeluLayer();

            lb = zeros(4, 4);
            ub = ones(4, 4);
            IS = ImageStar(lb, ub);

            % NN calls: reach(inSet, method, option, relaxFactor, dis_opt, lp_solver)
            IS_out = layer.reach(IS, 'approx-star', 'single', 0, [], 'linprog');

            testCase.verifyClass(IS_out, 'ImageStar');
        end

        function test_reach_nn_convention_6args(testCase)
            % Test reach with 6 arguments
            layer = GeluLayer();

            lb = zeros(4, 4);
            ub = ones(4, 4);
            IS = ImageStar(lb, ub);

            IS_out = layer.reach(IS, 'approx-star', 'single', 0, []);

            testCase.verifyClass(IS_out, 'ImageStar');
        end

        %% Soundness Tests

        function test_soundness_positive(testCase)
            % Verify output bounds contain actual GELU outputs
            layer = GeluLayer();

            lb = 0.5 * ones(10, 1);
            ub = 1.5 * ones(10, 1);

            [out_lb, out_ub] = layer.compute_gelu_bounds(lb, ub);

            % Sample random inputs and verify containment
            for i = 1:100
                x = lb + rand(10, 1) .* (ub - lb);
                y = layer.evaluate(x);

                testCase.verifyGreaterThanOrEqual(y, out_lb - 1e-6);
                testCase.verifyLessThanOrEqual(y, out_ub + 1e-6);
            end
        end

        function test_soundness_mixed(testCase)
            % Verify soundness with mixed bounds
            layer = GeluLayer();

            lb = -1 * ones(10, 1);
            ub = 1 * ones(10, 1);

            [out_lb, out_ub] = layer.compute_gelu_bounds(lb, ub);

            for i = 1:100
                x = lb + rand(10, 1) .* (ub - lb);
                y = layer.evaluate(x);

                testCase.verifyGreaterThanOrEqual(y, out_lb - 1e-6);
                testCase.verifyLessThanOrEqual(y, out_ub + 1e-6);
            end
        end

        %% Parse Tests

        function test_parse_matlab_layer(testCase)
            % Test parsing MATLAB geluLayer (if available)
            if ~exist('geluLayer', 'file')
                testCase.assumeFail('geluLayer not available in this MATLAB version');
            end

            matlab_layer = geluLayer('Name', 'test_gelu');
            layer = GeluLayer.parse(matlab_layer);

            testCase.verifyClass(layer, 'GeluLayer');
            testCase.verifyEqual(layer.Name, 'test_gelu');
        end

        %% Utility Tests

        function test_toGPU(testCase)
            % Test GPU conversion (should not error)
            layer = GeluLayer();
            layer_gpu = layer.toGPU();
            testCase.verifyClass(layer_gpu, 'GeluLayer');
        end

        function test_changeParamsPrecision(testCase)
            % Test precision change (should not error - no params)
            layer = GeluLayer();
            layer_single = layer.changeParamsPrecision('single');
            testCase.verifyClass(layer_single, 'GeluLayer');
        end

    end

end
