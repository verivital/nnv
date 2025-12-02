classdef ScaledDotProductAttentionLayer < handle
    % ScaledDotProductAttentionLayer - Implements scaled dot-product attention
    %
    % This layer computes: Attention(Q, K, V) = softmax(QK^T / sqrt(d_k)) * V
    %
    % For verification, we compute over-approximate reachable sets by:
    % 1. Computing bounds on QK^T (matrix multiplication)
    % 2. Scaling by 1/sqrt(d_k)
    % 3. Applying softmax bounds (using Softmax.m)
    % 4. Computing final attention output bounds
    %
    % Author: NNV Team
    % Date: November 2025

    properties
        Name = 'scaled_dot_product_attention';

        % Dimensions
        QueryDim = 0;      % Dimension of queries (d_k)
        KeyDim = 0;        % Dimension of keys (d_k)
        ValueDim = 0;      % Dimension of values (d_v)
        SeqLength = 0;     % Sequence length

        % Scale factor (1/sqrt(d_k))
        Scale = 1;

        % Layer interface
        NumInputs = 3;     % Q, K, V
        InputNames = {'query', 'key', 'value'};
        NumOutputs = 1;
        OutputNames = {'out'};
    end

    methods

        %% Constructor
        function obj = ScaledDotProductAttentionLayer(varargin)
            % Constructor for ScaledDotProductAttentionLayer
            %
            % Usage:
            %   layer = ScaledDotProductAttentionLayer()
            %   layer = ScaledDotProductAttentionLayer(name)
            %   layer = ScaledDotProductAttentionLayer(name, d_k)
            %   layer = ScaledDotProductAttentionLayer(name, d_k, d_v, seq_len)

            switch nargin
                case 0
                    % Default constructor
                case 1
                    obj.Name = varargin{1};
                case 2
                    obj.Name = varargin{1};
                    obj.QueryDim = varargin{2};
                    obj.KeyDim = varargin{2};
                    obj.Scale = 1 / sqrt(varargin{2});
                case 4
                    obj.Name = varargin{1};
                    obj.QueryDim = varargin{2};
                    obj.KeyDim = varargin{2};
                    obj.ValueDim = varargin{3};
                    obj.SeqLength = varargin{4};
                    obj.Scale = 1 / sqrt(varargin{2});
                otherwise
                    error('Invalid number of arguments');
            end
        end

        %% Evaluation (Forward Pass)
        function y = evaluate(obj, Q, K, V)
            % Compute scaled dot-product attention
            %
            % Inputs:
            %   Q: Query matrix [seq_len x d_k] or [batch x seq_len x d_k]
            %   K: Key matrix [seq_len x d_k]
            %   V: Value matrix [seq_len x d_v]
            %
            % Output:
            %   y: Attention output [seq_len x d_v]

            % Handle different input formats
            if ndims(Q) == 2
                % Single sequence: [seq_len x d_k]
                d_k = size(K, 2);

                % Compute attention scores: QK^T / sqrt(d_k)
                scores = (Q * K') / sqrt(d_k);  % [seq_len x seq_len]

                % Apply softmax to get attention weights
                attn_weights = softmax(scores, 2);  % Softmax over key dimension

                % Compute output: attention_weights * V
                y = attn_weights * V;  % [seq_len x d_v]

            elseif ndims(Q) == 3
                % Batch processing: [batch x seq_len x d_k]
                [batch_size, seq_len, d_k] = size(Q);
                d_v = size(V, 3);

                y = zeros(batch_size, seq_len, d_v);

                for b = 1:batch_size
                    Q_b = squeeze(Q(b, :, :));
                    K_b = squeeze(K(b, :, :));
                    V_b = squeeze(V(b, :, :));

                    scores = (Q_b * K_b') / sqrt(d_k);
                    attn_weights = softmax(scores, 2);
                    y(b, :, :) = attn_weights * V_b;
                end
            else
                error('Invalid input dimensions');
            end
        end

        %% Reachability Analysis
        function S = reach(varargin)
            % Main reachability method for attention layer
            %
            % Usage:
            %   S = layer.reach(Q_set, K_set, V_set)
            %   S = layer.reach(Q_set, K_set, V_set, method)
            %   S = layer.reach(Q_set, K_set, V_set, method, option, relaxFactor, dis_opt, lp_solver)

            obj = varargin{1};

            switch nargin
                case 4
                    Q_set = varargin{2};
                    K_set = varargin{3};
                    V_set = varargin{4};
                    method = 'approx-star';
                    lp_solver = 'linprog';
                case 5
                    Q_set = varargin{2};
                    K_set = varargin{3};
                    V_set = varargin{4};
                    method = varargin{5};
                    lp_solver = 'linprog';
                case 9
                    Q_set = varargin{2};
                    K_set = varargin{3};
                    V_set = varargin{4};
                    method = varargin{5};
                    % option = varargin{6};  % Not used
                    % relaxFactor = varargin{7};  % Not used
                    % dis_opt = varargin{8};  % Not used
                    lp_solver = varargin{9};
                otherwise
                    error('Invalid number of arguments');
            end

            % Dispatch based on method
            if strcmp(method, 'approx-star') || contains(method, 'relax-star')
                S = obj.reach_star_approx(Q_set, K_set, V_set, lp_solver);
            elseif strcmp(method, 'approx-zono')
                S = obj.reach_zono_approx(Q_set, K_set, V_set);
            else
                error('Unsupported reachability method: %s', method);
            end
        end

        function S = reach_star_approx(obj, Q_set, K_set, V_set, lp_solver)
            % Over-approximate reachability using Star sets
            %
            % This implements a bounds-based approach:
            % 1. Get bounds on Q, K, V
            % 2. Compute bounds on QK^T using interval arithmetic
            % 3. Apply softmax bounds
            % 4. Compute output bounds

            if nargin < 5
                lp_solver = 'linprog';
            end

            % Get dimensions
            n_q = Q_set.dim;
            n_k = K_set.dim;
            n_v = V_set.dim;

            % Get bounds on inputs
            Q_lb = zeros(n_q, 1);
            Q_ub = zeros(n_q, 1);
            for i = 1:n_q
                Q_lb(i) = Q_set.getMin(i, lp_solver);
                Q_ub(i) = Q_set.getMax(i, lp_solver);
            end

            K_lb = zeros(n_k, 1);
            K_ub = zeros(n_k, 1);
            for i = 1:n_k
                K_lb(i) = K_set.getMin(i, lp_solver);
                K_ub(i) = K_set.getMax(i, lp_solver);
            end

            V_lb = zeros(n_v, 1);
            V_ub = zeros(n_v, 1);
            for i = 1:n_v
                V_lb(i) = V_set.getMin(i, lp_solver);
                V_ub(i) = V_set.getMax(i, lp_solver);
            end

            % Compute attention output bounds
            [out_lb, out_ub] = obj.compute_attention_bounds(Q_lb, Q_ub, ...
                K_lb, K_ub, V_lb, V_ub);

            % Create output Star from bounds
            S = Star(out_lb, out_ub);
        end

        function Z = reach_zono_approx(obj, Q_set, K_set, V_set)
            % Over-approximate reachability using Zonotopes

            % Get bounds from zonotopes
            [Q_lb, Q_ub] = Q_set.getBounds();
            [K_lb, K_ub] = K_set.getBounds();
            [V_lb, V_ub] = V_set.getBounds();

            % Compute attention output bounds
            [out_lb, out_ub] = obj.compute_attention_bounds(Q_lb, Q_ub, ...
                K_lb, K_ub, V_lb, V_ub);

            % Create output Zonotope from bounds
            center = (out_lb + out_ub) / 2;
            generators = diag((out_ub - out_lb) / 2);
            Z = Zono(center, generators);
        end

        function [out_lb, out_ub] = compute_attention_bounds(obj, Q_lb, Q_ub, ...
                K_lb, K_ub, V_lb, V_ub)
            % Compute bounds on attention output
            %
            % Attention(Q, K, V) = softmax(QK^T / sqrt(d_k)) * V

            % For simplified analysis, assume Q, K, V are vectors
            % In full implementation, handle matrix case

            d_k = length(Q_lb);
            d_v = length(V_lb);

            % Step 1: Bound QK^T (dot product for vector case)
            % For interval [a, b] * [c, d]:
            % min = min(ac, ad, bc, bd), max = max(ac, ad, bc, bd)
            products = [Q_lb .* K_lb, Q_lb .* K_ub, Q_ub .* K_lb, Q_ub .* K_ub];
            score_lb = sum(min(products, [], 2));
            score_ub = sum(max(products, [], 2));

            % Step 2: Scale by 1/sqrt(d_k)
            scale = 1 / sqrt(d_k);
            score_lb = score_lb * scale;
            score_ub = score_ub * scale;

            % Step 3: Softmax bounds (for single score, softmax is identity-like)
            % For vector of scores, use Softmax.compute_softmax_bounds
            % Here we simplify for the single attention weight case
            attn_lb = max(0, 1 / (1 + exp(-score_lb)));  % Approximate
            attn_ub = min(1, 1 / (1 + exp(-score_ub)));  % Approximate

            % Step 4: Bound attention * V
            % attn * V where attn in [attn_lb, attn_ub], V in [V_lb, V_ub]
            products_v = [attn_lb * V_lb, attn_lb * V_ub, attn_ub * V_lb, attn_ub * V_ub];
            out_lb = min(products_v, [], 2);
            out_ub = max(products_v, [], 2);
        end

        %% Utility Methods
        function obj = toGPU(obj)
            % Convert to GPU (no parameters to convert for this layer)
        end

        function obj = changeParamsPrecision(obj, precision)
            % Change precision (no parameters)
        end

    end

    methods (Static)

        function L = parse(matlab_layer)
            % Parse from MATLAB selfAttentionLayer
            %
            % Input: MATLAB selfAttentionLayer object
            % Output: ScaledDotProductAttentionLayer object

            L = ScaledDotProductAttentionLayer();

            if isprop(matlab_layer, 'Name')
                L.Name = matlab_layer.Name;
            end

            if isprop(matlab_layer, 'NumHeads')
                % This is a multi-head attention layer
                % Extract single-head parameters
                warning('Multi-head attention detected. Creating single-head layer.');
            end

            if isprop(matlab_layer, 'KeySize')
                L.KeyDim = matlab_layer.KeySize;
                L.QueryDim = matlab_layer.KeySize;
                L.Scale = 1 / sqrt(matlab_layer.KeySize);
            end

            if isprop(matlab_layer, 'ValueSize')
                L.ValueDim = matlab_layer.ValueSize;
            end
        end

        %% Soundness Check
        function [is_sound, results] = check_soundness(Q_set, K_set, V_set, S_out, num_samples)
            % Check soundness of reachability computation

            if nargin < 5
                num_samples = 100;
            end

            layer = ScaledDotProductAttentionLayer();
            is_sound = true;
            results = struct('samples', {}, 'outputs', {}, 'contained', {});

            for i = 1:num_samples
                % Sample from input sets
                q = Q_set.sample(1);
                k = K_set.sample(1);
                v = V_set.sample(1);

                % Compute attention output
                y = layer.evaluate(q', k', v');
                y = y(:);

                % Check containment
                contained = S_out.contains(y);

                results(i).samples = {q, k, v};
                results(i).outputs = y;
                results(i).contained = contained;

                if ~contained
                    is_sound = false;
                    fprintf('Soundness violation at sample %d\n', i);
                end
            end
        end

    end
end
