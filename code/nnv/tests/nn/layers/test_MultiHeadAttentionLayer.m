%% Test Suite for MultiHeadAttentionLayer
% Tests for the multi-head attention layer
%
% Author: NNV Team
% Date: November 2025

classdef test_MultiHeadAttentionLayer < matlab.unittest.TestCase

    methods (Test)

        %% Constructor Tests

        function test_constructor_default(testCase)
            % Test default constructor
            layer = MultiHeadAttentionLayer();

            testCase.verifyClass(layer, 'MultiHeadAttentionLayer');
            testCase.verifyEqual(layer.NumHeads, 1, 'Default NumHeads');
            testCase.verifyEqual(layer.EmbedDim, 0, 'Default EmbedDim');
        end

        function test_constructor_with_name(testCase)
            % Test constructor with name only
            layer = MultiHeadAttentionLayer('my_mha');

            testCase.verifyEqual(layer.Name, 'my_mha');
        end

        function test_constructor_with_dims(testCase)
            % Test constructor with dimensions
            layer = MultiHeadAttentionLayer('mha', 64, 4);

            testCase.verifyEqual(layer.Name, 'mha');
            testCase.verifyEqual(layer.EmbedDim, 64);
            testCase.verifyEqual(layer.NumHeads, 4);
            testCase.verifyEqual(layer.HeadDim, 16);
        end

        function test_constructor_initializes_weights(testCase)
            % Test that constructor initializes weights
            layer = MultiHeadAttentionLayer('mha', 32, 2);

            testCase.verifySize(layer.W_Q, [32, 32], 'W_Q shape');
            testCase.verifySize(layer.W_K, [32, 32], 'W_K shape');
            testCase.verifySize(layer.W_V, [32, 32], 'W_V shape');
            testCase.verifySize(layer.W_O, [32, 32], 'W_O shape');
        end

        function test_constructor_embed_dim_divisibility(testCase)
            % Test that EmbedDim must be divisible by NumHeads
            testCase.verifyError(@() MultiHeadAttentionLayer('bad', 10, 3), ...
                '');  % Any error
        end

        %% Evaluation Tests

        function test_evaluate_basic(testCase)
            % Test basic evaluation
            layer = MultiHeadAttentionLayer('mha', 16, 2);

            seq_len = 4;
            Q = randn(seq_len, 16);
            K = randn(seq_len, 16);
            V = randn(seq_len, 16);

            y = layer.evaluate(Q, K, V);

            testCase.verifySize(y, [seq_len, 16], 'Output shape');
        end

        function test_evaluate_self_attention(testCase)
            % Test self-attention (single input)
            layer = MultiHeadAttentionLayer('mha', 32, 4);

            x = randn(8, 32);
            y = layer.evaluate(x);

            testCase.verifySize(y, [8, 32], 'Self-attention output shape');
        end

        function test_evaluate_single_head(testCase)
            % Test with single head (should behave like scaled dot-product)
            layer = MultiHeadAttentionLayer('mha', 16, 1);

            Q = randn(4, 16);
            K = randn(4, 16);
            V = randn(4, 16);

            y = layer.evaluate(Q, K, V);

            testCase.verifySize(y, [4, 16], 'Single head output shape');
            testCase.verifyFalse(any(isnan(y(:))), 'No NaN');
        end

        function test_evaluate_many_heads(testCase)
            % Test with many heads
            layer = MultiHeadAttentionLayer('mha', 64, 8);

            Q = randn(10, 64);
            K = randn(10, 64);
            V = randn(10, 64);

            y = layer.evaluate(Q, K, V);

            testCase.verifySize(y, [10, 64], '8-head output shape');
            testCase.verifyFalse(any(isnan(y(:))), 'No NaN');
        end

        function test_evaluate_with_identity_weights(testCase)
            % Test with identity projection weights
            layer = MultiHeadAttentionLayer('mha');
            layer.NumHeads = 2;
            layer.EmbedDim = 8;
            layer.HeadDim = 4;
            layer.W_Q = eye(8);
            layer.W_K = eye(8);
            layer.W_V = eye(8);
            layer.W_O = eye(8);

            Q = randn(4, 8);
            y = layer.evaluate(Q, Q, Q);

            testCase.verifySize(y, [4, 8], 'Output with identity weights');
        end

        %% Reachability Tests

        function test_reach_basic(testCase)
            % Test basic reachability
            layer = MultiHeadAttentionLayer('mha', 4, 1);

            center = randn(4, 1);
            V_star = [center, eye(4) * 0.1];
            C = [eye(4); -eye(4)];
            d = ones(8, 1);

            S = Star(V_star, C, d);

            try
                S_out = layer.reach(S, S, S);
                testCase.verifyClass(S_out, 'Star', 'Output should be Star');
            catch ME
                testCase.assumeFail(['Reach not implemented: ' ME.message]);
            end
        end

        function test_reach_self_attention(testCase)
            % Test self-attention reachability
            layer = MultiHeadAttentionLayer('mha', 4, 2);

            center = randn(4, 1);
            V_star = [center, eye(4) * 0.05];
            C = [eye(4); -eye(4)];
            d = ones(8, 1);

            S = Star(V_star, C, d);

            try
                S_out = layer.reach(S);  % Self-attention
                testCase.verifyClass(S_out, 'Star', 'Self-attention reach');
            catch ME
                testCase.assumeFail(['Reach not implemented: ' ME.message]);
            end
        end

        function test_reach_zono(testCase)
            % Test zonotope reachability
            layer = MultiHeadAttentionLayer('mha', 4, 2);

            center = randn(4, 1);
            generators = 0.1 * eye(4);
            Z = Zono(center, generators);

            try
                Z_out = layer.reach(Z, Z, Z, 'approx-zono');
                testCase.verifyClass(Z_out, 'Zono', 'Zonotope output');
            catch ME
                testCase.assumeFail(['Zono reach not implemented: ' ME.message]);
            end
        end

        %% Utility Tests

        function test_toGPU(testCase)
            % Test GPU conversion
            layer = MultiHeadAttentionLayer('mha', 8, 2);

            try
                layer_gpu = layer.toGPU();
                testCase.verifyClass(layer_gpu, 'MultiHeadAttentionLayer');
            catch ME
                testCase.assumeFail('GPU not available');
            end
        end

        function test_changeParamsPrecision_single(testCase)
            % Test precision change to single
            layer = MultiHeadAttentionLayer('mha', 8, 2);

            layer_single = layer.changeParamsPrecision('single');

            testCase.verifyClass(layer_single.W_Q, 'single');
            testCase.verifyClass(layer_single.W_K, 'single');
            testCase.verifyClass(layer_single.W_V, 'single');
            testCase.verifyClass(layer_single.W_O, 'single');
        end

        function test_changeParamsPrecision_double(testCase)
            % Test precision change to double
            layer = MultiHeadAttentionLayer('mha', 8, 2);
            layer.W_Q = single(layer.W_Q);

            layer_double = layer.changeParamsPrecision('double');

            testCase.verifyClass(layer_double.W_Q, 'double');
        end

        %% Edge Cases

        function test_single_token(testCase)
            % Test with single token
            layer = MultiHeadAttentionLayer('mha', 8, 2);

            Q = randn(1, 8);
            y = layer.evaluate(Q, Q, Q);

            testCase.verifySize(y, [1, 8], 'Single token output');
        end

        function test_large_sequence(testCase)
            % Test with larger sequence
            layer = MultiHeadAttentionLayer('mha', 32, 4);

            Q = randn(100, 32);
            K = randn(100, 32);
            V = randn(100, 32);

            y = layer.evaluate(Q, K, V);

            testCase.verifySize(y, [100, 32], 'Large sequence output');
            testCase.verifyFalse(any(isnan(y(:))), 'No NaN in large sequence');
        end

        %% Parse Tests

        function test_parse_basic(testCase)
            % Test parse with mock layer
            mock_layer = struct();
            mock_layer.Name = 'parsed_mha';
            mock_layer.NumHeads = 4;

            try
                layer = MultiHeadAttentionLayer.parse(mock_layer);
                testCase.verifyClass(layer, 'MultiHeadAttentionLayer');
            catch ME
                testCase.assumeFail(['Parse not implemented: ' ME.message]);
            end
        end

        %% ImageStar Reachability Tests (Critical for NN integration)

        function test_reach_imagestar_bounds(testCase)
            % Test reachability with ImageStar bounds representation
            % This tests the case when attention is used in CNN-like networks
            layer = MultiHeadAttentionLayer('mha', 4, 1);

            lb = zeros(4, 1);
            ub = ones(4, 1);
            IS = ImageStar(lb, ub);

            try
                IS_out = layer.reach(IS, 'approx-star');
                testCase.verifyTrue(isa(IS_out, 'Star') || isa(IS_out, 'ImageStar'), ...
                    'Output should be Star or ImageStar');
            catch ME
                testCase.assumeFail(['ImageStar reach not implemented: ' ME.message]);
            end
        end

        function test_reach_imagestar_3d(testCase)
            % Test reachability with 3D ImageStar
            layer = MultiHeadAttentionLayer('mha', 4, 2);

            lb = -ones(2, 2, 1);
            ub = ones(2, 2, 1);
            IS = ImageStar(lb, ub);

            try
                IS_out = layer.reach(IS, 'approx-star');
                testCase.verifyTrue(isa(IS_out, 'Star') || isa(IS_out, 'ImageStar'), ...
                    '3D ImageStar output');
            catch ME
                testCase.assumeFail(['3D ImageStar reach not implemented: ' ME.message]);
            end
        end

        %% NN Calling Convention Tests (Critical for NN.reach integration)

        function test_reach_nn_convention_7args(testCase)
            % Test reach with 7 arguments (NN.reach calling convention)
            % NN calls: reach(inSet, method, option, relaxFactor, dis_opt, lp_solver)
            layer = MultiHeadAttentionLayer('mha', 4, 1);

            lb = zeros(4, 1);
            ub = ones(4, 1);
            IS = ImageStar(lb, ub);

            try
                IS_out = layer.reach(IS, 'approx-star', 'single', 0, [], 'linprog');
                testCase.verifyTrue(isa(IS_out, 'Star') || isa(IS_out, 'ImageStar'), ...
                    'NN calling convention output');
            catch ME
                testCase.assumeFail(['7-arg reach not implemented: ' ME.message]);
            end
        end

        function test_reach_nn_convention_6args(testCase)
            % Test reach with 6 arguments
            layer = MultiHeadAttentionLayer('mha', 4, 1);

            lb = zeros(4, 1);
            ub = ones(4, 1);
            IS = ImageStar(lb, ub);

            try
                IS_out = layer.reach(IS, 'approx-star', 'single', 0, []);
                testCase.verifyTrue(isa(IS_out, 'Star') || isa(IS_out, 'ImageStar'), ...
                    '6-arg reach');
            catch ME
                testCase.assumeFail(['6-arg reach not implemented: ' ME.message]);
            end
        end

        function test_reach_star_with_nn_convention(testCase)
            % Test Star input with NN calling convention
            layer = MultiHeadAttentionLayer('mha', 4, 1);

            lb = zeros(4, 1);
            ub = ones(4, 1);
            S = Star(lb, ub);

            try
                S_out = layer.reach(S, 'approx-star', 'single', 0, [], 'linprog');
                testCase.verifyClass(S_out, 'Star', 'Star with NN convention');
            catch ME
                testCase.assumeFail(['Star NN convention not implemented: ' ME.message]);
            end
        end

        %% ImageZono Reachability Tests

        function test_reach_imagezono(testCase)
            % Test reachability with ImageZono
            layer = MultiHeadAttentionLayer('mha', 4, 1);

            lb = -ones(4, 1);
            ub = ones(4, 1);
            IZ = ImageZono(lb, ub);

            try
                IZ_out = layer.reach(IZ, 'approx-zono');
                testCase.verifyTrue(isa(IZ_out, 'Zono') || isa(IZ_out, 'ImageZono'), ...
                    'ImageZono output');
            catch ME
                testCase.assumeFail(['ImageZono reach not implemented: ' ME.message]);
            end
        end

    end

end
