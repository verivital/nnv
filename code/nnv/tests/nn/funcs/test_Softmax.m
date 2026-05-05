%% Test Suite for Softmax Reachability
% Tests for the Softmax class implementing softmax reachability analysis
%
% Author: NNV Team
% Date: November 2025

classdef test_Softmax < matlab.unittest.TestCase

    methods (Test)

        %% Basic Evaluation Tests

        function test_evaluate_basic(testCase)
            % Test basic softmax evaluation
            x = [1; 2; 3];
            y = Softmax.evaluate(x);

            % Check output is valid probability distribution
            testCase.verifyEqual(sum(y), 1, 'AbsTol', 1e-10, 'Softmax should sum to 1');
            testCase.verifyTrue(all(y >= 0), 'Softmax outputs should be non-negative');
            testCase.verifyTrue(all(y <= 1), 'Softmax outputs should be <= 1');
        end

        function test_evaluate_numerical_stability(testCase)
            % Test numerical stability with large values
            x = [1000; 1001; 1002];
            y = Softmax.evaluate(x);

            testCase.verifyEqual(sum(y), 1, 'AbsTol', 1e-10, 'Should handle large values');
            testCase.verifyFalse(any(isnan(y)), 'Should not produce NaN');
            testCase.verifyFalse(any(isinf(y)), 'Should not produce Inf');
        end

        function test_evaluate_uniform(testCase)
            % Test uniform input produces uniform output
            x = [5; 5; 5; 5];
            y = Softmax.evaluate(x);

            expected = ones(4, 1) / 4;
            testCase.verifyEqual(y, expected, 'AbsTol', 1e-10, 'Equal inputs should give uniform output');
        end

        function test_evaluate_ordering(testCase)
            % Test that larger input gives larger probability
            x = [1; 2; 3];
            y = Softmax.evaluate(x);

            testCase.verifyTrue(y(3) > y(2), 'Larger input should have larger probability');
            testCase.verifyTrue(y(2) > y(1), 'Larger input should have larger probability');
        end

        %% Bound Computation Tests

        function test_compute_bounds_basic(testCase)
            % Test basic bound computation
            lb = [0; 0; 0];
            ub = [1; 1; 1];

            [sm_lb, sm_ub] = Softmax.compute_softmax_bounds(lb, ub);

            % Check bounds are valid
            testCase.verifyTrue(all(sm_lb >= 0), 'Lower bounds should be >= 0');
            testCase.verifyTrue(all(sm_ub <= 1), 'Upper bounds should be <= 1');
            testCase.verifyTrue(all(sm_lb <= sm_ub), 'Lower bounds should be <= upper bounds');
        end

        function test_compute_bounds_soundness(testCase)
            % Test soundness: all sampled points should have softmax within bounds
            lb = [-1; 0; 1];
            ub = [0; 1; 2];

            [sm_lb, sm_ub] = Softmax.compute_softmax_bounds(lb, ub);

            % Sample many points and verify containment
            num_samples = 100;
            for i = 1:num_samples
                x = lb + (ub - lb) .* rand(3, 1);
                y = Softmax.evaluate(x);

                testCase.verifyTrue(all(y >= sm_lb - 1e-10), ...
                    sprintf('Sample %d softmax should be >= lower bounds', i));
                testCase.verifyTrue(all(y <= sm_ub + 1e-10), ...
                    sprintf('Sample %d softmax should be <= upper bounds', i));
            end
        end

        function test_compute_bounds_tight_input(testCase)
            % Test with tight (point) input
            x = [1; 2; 3];
            [sm_lb, sm_ub] = Softmax.compute_softmax_bounds(x, x);

            y = Softmax.evaluate(x);

            % Bounds should be close to actual softmax
            testCase.verifyEqual(sm_lb, y, 'AbsTol', 0.01, 'Tight bounds should approximate softmax');
            testCase.verifyEqual(sm_ub, y, 'AbsTol', 0.01, 'Tight bounds should approximate softmax');
        end

        function test_compute_bounds_wide_input(testCase)
            % Test with wide input range
            lb = [-5; -5; -5];
            ub = [5; 5; 5];

            [sm_lb, sm_ub] = Softmax.compute_softmax_bounds(lb, ub);

            % Bounds should be valid but potentially loose
            testCase.verifyTrue(all(sm_lb >= 0), 'Lower bounds should be >= 0');
            testCase.verifyTrue(all(sm_ub <= 1), 'Upper bounds should be <= 1');
        end

        %% Star Set Reachability Tests

        function test_reach_star_approx_basic(testCase)
            % Test basic star set reachability
            center = [1; 2; 3];
            V = [center, eye(3) * 0.1];
            C = [eye(3); -eye(3)];
            d = ones(6, 1);

            S_in = Star(V, C, d);
            S_out = Softmax.reach_star_approx(S_in);

            testCase.verifyClass(S_out, 'Star', 'Output should be a Star set');
            testCase.verifyEqual(S_out.dim, 3, 'Output dimension should match input');
        end

        function test_reach_star_approx_soundness(testCase)
            % Test soundness of star reachability
            center = [0; 0; 0];
            V = [center, eye(3) * 0.5];
            C = [eye(3); -eye(3)];
            d = ones(6, 1);

            S_in = Star(V, C, d);
            S_out = Softmax.reach_star_approx(S_in);

            % Sample from input and verify output containment
            num_samples = 50;
            for i = 1:num_samples
                % Sample from input star
                alpha = 2 * rand(3, 1) - 1;  % Random in [-1, 1]^3
                x = center + V(:, 2:end) * alpha;

                % Compute softmax
                y = Softmax.evaluate(x);

                % Check containment (approximately)
                [lb, ub] = S_out.getRanges();
                testCase.verifyTrue(all(y >= lb - 0.1), ...
                    sprintf('Sample %d should be within output bounds', i));
                testCase.verifyTrue(all(y <= ub + 0.1), ...
                    sprintf('Sample %d should be within output bounds', i));
            end
        end

        %% Zonotope Reachability Tests

        function test_reach_zono_approx_basic(testCase)
            % Test basic zonotope reachability
            center = [1; 2; 3];
            generators = eye(3) * 0.1;

            Z_in = Zono(center, generators);
            Z_out = Softmax.reach_zono_approx(Z_in);

            testCase.verifyClass(Z_out, 'Zono', 'Output should be a Zonotope');
        end

        %% Edge Cases

        function test_single_element(testCase)
            % Test single element input
            x = 5;
            y = Softmax.evaluate(x);

            testCase.verifyEqual(y, 1, 'AbsTol', 1e-10, 'Single element softmax should be 1');
        end

        function test_two_elements(testCase)
            % Test two element input (binary classification)
            x = [0; 1];
            y = Softmax.evaluate(x);

            testCase.verifyEqual(sum(y), 1, 'AbsTol', 1e-10, 'Should sum to 1');
            testCase.verifyTrue(y(2) > y(1), 'Larger input should have larger probability');
        end

        function test_negative_inputs(testCase)
            % Test negative inputs
            x = [-3; -2; -1];
            y = Softmax.evaluate(x);

            testCase.verifyEqual(sum(y), 1, 'AbsTol', 1e-10, 'Should handle negative inputs');
            testCase.verifyTrue(all(y > 0), 'All outputs should be positive');
        end

    end

end
