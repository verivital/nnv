classdef MultiHeadAttentionLayer < handle
    % MultiHeadAttentionLayer - Implements multi-head attention mechanism
    %
    % This layer computes multi-head attention as in "Attention Is All You Need":
    %   MultiHead(Q, K, V) = Concat(head_1, ..., head_h) * W_O
    %   where head_i = Attention(Q * W_Q^i, K * W_K^i, V * W_V^i)
    %
    % For verification, we compute over-approximate reachable sets by:
    % 1. Computing bounds through linear projections for each head
    % 2. Computing attention bounds for each head
    % 3. Concatenating and projecting output bounds
    %
    % Author: NNV Team
    % Date: November 2025

    properties
        Name = 'multi_head_attention';

        % Dimensions
        NumHeads = 1;          % Number of attention heads
        EmbedDim = 0;          % Input embedding dimension (d_model)
        HeadDim = 0;           % Dimension per head (d_k = d_model / num_heads)

        % Projection weights (learned)
        W_Q = [];              % Query projection [EmbedDim x EmbedDim]
        W_K = [];              % Key projection [EmbedDim x EmbedDim]
        W_V = [];              % Value projection [EmbedDim x EmbedDim]
        W_O = [];              % Output projection [EmbedDim x EmbedDim]

        % Biases (optional)
        b_Q = [];
        b_K = [];
        b_V = [];
        b_O = [];

        % Layer interface
        NumInputs = 3;         % Q, K, V (or 1 for self-attention)
        InputNames = {'query', 'key', 'value'};
        NumOutputs = 1;
        OutputNames = {'out'};
    end

    methods

        %% Constructor
        function obj = MultiHeadAttentionLayer(varargin)
            % Constructor for MultiHeadAttentionLayer
            %
            % Usage:
            %   layer = MultiHeadAttentionLayer()
            %   layer = MultiHeadAttentionLayer(name)
            %   layer = MultiHeadAttentionLayer(name, embed_dim, num_heads)
            %   layer = MultiHeadAttentionLayer(name, embed_dim, num_heads, W_Q, W_K, W_V, W_O)

            switch nargin
                case 0
                    % Default constructor
                case 1
                    obj.Name = varargin{1};
                case 3
                    obj.Name = varargin{1};
                    obj.EmbedDim = varargin{2};
                    obj.NumHeads = varargin{3};
                    obj.HeadDim = obj.EmbedDim / obj.NumHeads;

                    if mod(obj.EmbedDim, obj.NumHeads) ~= 0
                        error('EmbedDim must be divisible by NumHeads');
                    end

                    % Initialize random weights (Xavier initialization)
                    scale = sqrt(2 / (obj.EmbedDim + obj.HeadDim));
                    obj.W_Q = randn(obj.EmbedDim, obj.EmbedDim) * scale;
                    obj.W_K = randn(obj.EmbedDim, obj.EmbedDim) * scale;
                    obj.W_V = randn(obj.EmbedDim, obj.EmbedDim) * scale;
                    obj.W_O = randn(obj.EmbedDim, obj.EmbedDim) * scale;

                case 7
                    obj.Name = varargin{1};
                    obj.EmbedDim = varargin{2};
                    obj.NumHeads = varargin{3};
                    obj.HeadDim = obj.EmbedDim / obj.NumHeads;
                    obj.W_Q = varargin{4};
                    obj.W_K = varargin{5};
                    obj.W_V = varargin{6};
                    obj.W_O = varargin{7};

                otherwise
                    error('Invalid number of arguments');
            end
        end

        %% Evaluation (Forward Pass)
        function y = evaluate(obj, Q, K, V)
            % Compute multi-head attention
            %
            % Inputs:
            %   Q: Query matrix [seq_len x embed_dim]
            %   K: Key matrix [seq_len x embed_dim]
            %   V: Value matrix [seq_len x embed_dim]
            %
            % Output:
            %   y: Attention output [seq_len x embed_dim]

            if nargin == 2
                % Self-attention: Q = K = V
                K = Q;
                V = Q;
            end

            [seq_len, embed_dim] = size(Q);

            if embed_dim ~= obj.EmbedDim && obj.EmbedDim > 0
                error('Input dimension mismatch');
            end

            % If weights not initialized, use identity (pass-through)
            if isempty(obj.W_Q)
                obj.EmbedDim = embed_dim;
                obj.HeadDim = embed_dim / obj.NumHeads;
                obj.W_Q = eye(embed_dim);
                obj.W_K = eye(embed_dim);
                obj.W_V = eye(embed_dim);
                obj.W_O = eye(embed_dim);
            end

            % Project Q, K, V
            Q_proj = Q * obj.W_Q;  % [seq_len x embed_dim]
            K_proj = K * obj.W_K;
            V_proj = V * obj.W_V;

            % Add biases if present
            if ~isempty(obj.b_Q)
                Q_proj = Q_proj + obj.b_Q';
            end
            if ~isempty(obj.b_K)
                K_proj = K_proj + obj.b_K';
            end
            if ~isempty(obj.b_V)
                V_proj = V_proj + obj.b_V';
            end

            % Reshape for multi-head: [seq_len x num_heads x head_dim]
            Q_heads = reshape(Q_proj, [seq_len, obj.NumHeads, obj.HeadDim]);
            K_heads = reshape(K_proj, [seq_len, obj.NumHeads, obj.HeadDim]);
            V_heads = reshape(V_proj, [seq_len, obj.NumHeads, obj.HeadDim]);

            % Compute attention for each head
            head_outputs = zeros(seq_len, obj.NumHeads, obj.HeadDim);

            for h = 1:obj.NumHeads
                Q_h = squeeze(Q_heads(:, h, :));  % [seq_len x head_dim]
                K_h = squeeze(K_heads(:, h, :));
                V_h = squeeze(V_heads(:, h, :));

                % Scaled dot-product attention
                d_k = obj.HeadDim;
                scores = (Q_h * K_h') / sqrt(d_k);  % [seq_len x seq_len]
                attn_weights = softmax(scores, 2);  % Softmax over keys
                head_outputs(:, h, :) = attn_weights * V_h;
            end

            % Concatenate heads: [seq_len x embed_dim]
            concat_output = reshape(head_outputs, [seq_len, obj.EmbedDim]);

            % Output projection
            y = concat_output * obj.W_O;

            if ~isempty(obj.b_O)
                y = y + obj.b_O';
            end
        end

        %% Reachability Analysis
        function S = reach(varargin)
            % Main reachability method for multi-head attention
            %
            % Usage:
            %   S = layer.reach(input_set)           % Self-attention
            %   S = layer.reach(Q_set, K_set, V_set) % Cross-attention
            %   S = layer.reach(..., method)
            %   S = layer.reach(..., method, option, relaxFactor, dis_opt, lp_solver)

            obj = varargin{1};

            switch nargin
                case 2
                    % Self-attention
                    Q_set = varargin{2};
                    K_set = Q_set;
                    V_set = Q_set;
                    method = 'approx-star';
                    lp_solver = 'linprog';
                case 3
                    Q_set = varargin{2};
                    K_set = Q_set;
                    V_set = Q_set;
                    method = varargin{3};
                    lp_solver = 'linprog';
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
            % 1. Get bounds on projected Q, K, V
            % 2. Compute attention bounds for each head
            % 3. Compute output bounds

            if nargin < 5
                lp_solver = 'linprog';
            end

            % Get bounds on inputs
            n = Q_set.dim;
            Q_lb = zeros(n, 1);
            Q_ub = zeros(n, 1);
            K_lb = zeros(n, 1);
            K_ub = zeros(n, 1);
            V_lb = zeros(n, 1);
            V_ub = zeros(n, 1);

            for i = 1:n
                Q_lb(i) = Q_set.getMin(i, lp_solver);
                Q_ub(i) = Q_set.getMax(i, lp_solver);
                K_lb(i) = K_set.getMin(i, lp_solver);
                K_ub(i) = K_set.getMax(i, lp_solver);
                V_lb(i) = V_set.getMin(i, lp_solver);
                V_ub(i) = V_set.getMax(i, lp_solver);
            end

            % Compute output bounds through multi-head attention
            [out_lb, out_ub] = obj.compute_mha_bounds(Q_lb, Q_ub, K_lb, K_ub, V_lb, V_ub);

            % Create output Star from bounds
            S = Star(out_lb, out_ub);
        end

        function Z = reach_zono_approx(obj, Q_set, K_set, V_set)
            % Over-approximate reachability using Zonotopes

            % Get bounds from zonotopes
            [Q_lb, Q_ub] = Q_set.getBounds();
            [K_lb, K_ub] = K_set.getBounds();
            [V_lb, V_ub] = V_set.getBounds();

            % Compute output bounds
            [out_lb, out_ub] = obj.compute_mha_bounds(Q_lb, Q_ub, K_lb, K_ub, V_lb, V_ub);

            % Create output Zonotope from bounds
            center = (out_lb + out_ub) / 2;
            generators = diag((out_ub - out_lb) / 2);
            Z = Zono(center, generators);
        end

        function [out_lb, out_ub] = compute_mha_bounds(obj, Q_lb, Q_ub, K_lb, K_ub, V_lb, V_ub)
            % Compute bounds on multi-head attention output
            %
            % Uses interval arithmetic through:
            % 1. Linear projections
            % 2. Per-head attention
            % 3. Output projection

            n = length(Q_lb);

            % If weights not set, use identity
            if isempty(obj.W_Q)
                obj.W_Q = eye(n);
                obj.W_K = eye(n);
                obj.W_V = eye(n);
                obj.W_O = eye(n);
                obj.EmbedDim = n;
                obj.HeadDim = n / obj.NumHeads;
            end

            % Project inputs through linear layers (interval arithmetic)
            [Q_proj_lb, Q_proj_ub] = obj.interval_matmul(Q_lb, Q_ub, obj.W_Q);
            [K_proj_lb, K_proj_ub] = obj.interval_matmul(K_lb, K_ub, obj.W_K);
            [V_proj_lb, V_proj_ub] = obj.interval_matmul(V_lb, V_ub, obj.W_V);

            % Compute attention bounds for each head
            head_lb = zeros(n, 1);
            head_ub = zeros(n, 1);

            head_size = obj.HeadDim;
            for h = 1:obj.NumHeads
                idx_start = (h-1) * head_size + 1;
                idx_end = h * head_size;

                q_h_lb = Q_proj_lb(idx_start:idx_end);
                q_h_ub = Q_proj_ub(idx_start:idx_end);
                k_h_lb = K_proj_lb(idx_start:idx_end);
                k_h_ub = K_proj_ub(idx_start:idx_end);
                v_h_lb = V_proj_lb(idx_start:idx_end);
                v_h_ub = V_proj_ub(idx_start:idx_end);

                [h_lb, h_ub] = obj.compute_attention_bounds_single_head(...
                    q_h_lb, q_h_ub, k_h_lb, k_h_ub, v_h_lb, v_h_ub);

                head_lb(idx_start:idx_end) = h_lb;
                head_ub(idx_start:idx_end) = h_ub;
            end

            % Output projection
            [out_lb, out_ub] = obj.interval_matmul(head_lb, head_ub, obj.W_O);
        end

        function [y_lb, y_ub] = interval_matmul(~, x_lb, x_ub, W)
            % Interval matrix multiplication: y = x * W
            % For y_i = sum_j x_j * W_ji

            n_out = size(W, 2);
            y_lb = zeros(n_out, 1);
            y_ub = zeros(n_out, 1);

            for i = 1:n_out
                w_col = W(:, i);

                % For each output, compute bounds of dot product
                products = [x_lb .* w_col, x_lb .* w_col, x_ub .* w_col, x_ub .* w_col];

                % Handle sign cases properly
                pos_mask = w_col >= 0;
                neg_mask = ~pos_mask;

                y_lb(i) = sum(x_lb(pos_mask) .* w_col(pos_mask)) + ...
                          sum(x_ub(neg_mask) .* w_col(neg_mask));
                y_ub(i) = sum(x_ub(pos_mask) .* w_col(pos_mask)) + ...
                          sum(x_lb(neg_mask) .* w_col(neg_mask));
            end
        end

        function [out_lb, out_ub] = compute_attention_bounds_single_head(obj, ...
                q_lb, q_ub, k_lb, k_ub, v_lb, v_ub)
            % Compute bounds for a single attention head

            d_k = length(q_lb);

            % Bound QK^T (dot product)
            products = [q_lb .* k_lb, q_lb .* k_ub, q_ub .* k_lb, q_ub .* k_ub];
            score_lb = sum(min(products, [], 2));
            score_ub = sum(max(products, [], 2));

            % Scale by 1/sqrt(d_k)
            scale = 1 / sqrt(d_k);
            score_lb = score_lb * scale;
            score_ub = score_ub * scale;

            % Softmax bounds (approximation for single score)
            attn_lb = 1 / (1 + exp(-score_lb));
            attn_ub = 1 / (1 + exp(-score_ub));
            attn_lb = max(0, min(1, attn_lb));
            attn_ub = max(0, min(1, attn_ub));

            % Bound attention * V
            products_v = [attn_lb * v_lb, attn_lb * v_ub, attn_ub * v_lb, attn_ub * v_ub];
            out_lb = min(products_v, [], 2);
            out_ub = max(products_v, [], 2);
        end

        %% Utility Methods
        function obj = toGPU(obj)
            % Convert weights to GPU arrays
            if ~isempty(obj.W_Q)
                obj.W_Q = gpuArray(obj.W_Q);
                obj.W_K = gpuArray(obj.W_K);
                obj.W_V = gpuArray(obj.W_V);
                obj.W_O = gpuArray(obj.W_O);
            end
            if ~isempty(obj.b_Q)
                obj.b_Q = gpuArray(obj.b_Q);
                obj.b_K = gpuArray(obj.b_K);
                obj.b_V = gpuArray(obj.b_V);
                obj.b_O = gpuArray(obj.b_O);
            end
        end

        function obj = changeParamsPrecision(obj, precision)
            % Change precision of parameters
            if ~isempty(obj.W_Q)
                if strcmp(precision, 'single')
                    obj.W_Q = single(obj.W_Q);
                    obj.W_K = single(obj.W_K);
                    obj.W_V = single(obj.W_V);
                    obj.W_O = single(obj.W_O);
                elseif strcmp(precision, 'double')
                    obj.W_Q = double(obj.W_Q);
                    obj.W_K = double(obj.W_K);
                    obj.W_V = double(obj.W_V);
                    obj.W_O = double(obj.W_O);
                end
            end
            if ~isempty(obj.b_Q)
                if strcmp(precision, 'single')
                    obj.b_Q = single(obj.b_Q);
                    obj.b_K = single(obj.b_K);
                    obj.b_V = single(obj.b_V);
                    obj.b_O = single(obj.b_O);
                elseif strcmp(precision, 'double')
                    obj.b_Q = double(obj.b_Q);
                    obj.b_K = double(obj.b_K);
                    obj.b_V = double(obj.b_V);
                    obj.b_O = double(obj.b_O);
                end
            end
        end

    end

    methods (Static)

        function L = parse(matlab_layer)
            % Parse from MATLAB selfAttentionLayer
            %
            % Input: MATLAB selfAttentionLayer object
            % Output: MultiHeadAttentionLayer object
            %
            % MATLAB selfAttentionLayer properties:
            %   - NumHeads: Number of attention heads
            %   - NumKeyValueChannels: Channel dimension for K and V
            %   - NumQueryChannels: Channel dimension for Q (optional)
            %   - QueryWeights, KeyWeights, ValueWeights, OutputWeights
            %   - QueryBias, KeyBias, ValueBias, OutputBias (if HasBias)

            L = MultiHeadAttentionLayer();

            % Parse name
            if isprop(matlab_layer, 'Name')
                L.Name = matlab_layer.Name;
            end

            % Parse number of heads
            if isprop(matlab_layer, 'NumHeads')
                L.NumHeads = matlab_layer.NumHeads;
            end

            % Parse dimensions
            if isprop(matlab_layer, 'NumChannels')
                L.EmbedDim = matlab_layer.NumChannels;
            elseif isprop(matlab_layer, 'NumKeyValueChannels')
                L.EmbedDim = matlab_layer.NumKeyValueChannels;
            end

            if L.NumHeads > 0 && L.EmbedDim > 0
                L.HeadDim = L.EmbedDim / L.NumHeads;
            end

            % Extract projection weights if available
            % MATLAB stores weights as [OutputSize x InputSize]
            if isprop(matlab_layer, 'QueryWeights') && ~isempty(matlab_layer.QueryWeights)
                L.W_Q = double(matlab_layer.QueryWeights)';  % Transpose for our convention
            end
            if isprop(matlab_layer, 'KeyWeights') && ~isempty(matlab_layer.KeyWeights)
                L.W_K = double(matlab_layer.KeyWeights)';
            end
            if isprop(matlab_layer, 'ValueWeights') && ~isempty(matlab_layer.ValueWeights)
                L.W_V = double(matlab_layer.ValueWeights)';
            end
            if isprop(matlab_layer, 'OutputWeights') && ~isempty(matlab_layer.OutputWeights)
                L.W_O = double(matlab_layer.OutputWeights)';
            end

            % Extract biases if available
            if isprop(matlab_layer, 'QueryBias') && ~isempty(matlab_layer.QueryBias)
                L.b_Q = double(matlab_layer.QueryBias);
            end
            if isprop(matlab_layer, 'KeyBias') && ~isempty(matlab_layer.KeyBias)
                L.b_K = double(matlab_layer.KeyBias);
            end
            if isprop(matlab_layer, 'ValueBias') && ~isempty(matlab_layer.ValueBias)
                L.b_V = double(matlab_layer.ValueBias);
            end
            if isprop(matlab_layer, 'OutputBias') && ~isempty(matlab_layer.OutputBias)
                L.b_O = double(matlab_layer.OutputBias);
            end

            % Fallback: try generic Weights property (older MATLAB versions)
            if isempty(L.W_Q) && isprop(matlab_layer, 'Weights')
                w = matlab_layer.Weights;
                if ~isempty(w)
                    d = L.EmbedDim;
                    if d > 0 && numel(w) >= 4*d*d
                        L.W_Q = reshape(w(1:d*d), [d, d]);
                        L.W_K = reshape(w(d*d+1:2*d*d), [d, d]);
                        L.W_V = reshape(w(2*d*d+1:3*d*d), [d, d]);
                        L.W_O = reshape(w(3*d*d+1:4*d*d), [d, d]);
                    end
                end
            end

            % Warning if weights not found (layer will use identity)
            if isempty(L.W_Q)
                warning('NNV:MultiHeadAttentionLayer:NoWeights', ...
                    'Could not extract attention weights. Using identity projections.');
            end
        end

        %% Soundness Check
        function [is_sound, results] = check_soundness(Q_set, K_set, V_set, S_out, num_samples)
            % Check soundness of reachability computation

            if nargin < 5
                num_samples = 100;
            end

            layer = MultiHeadAttentionLayer();
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
