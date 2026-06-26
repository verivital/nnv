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
                case {6, 7, 8}
                    % [46] Accept the 6/7/8-arg shapes (Q,K,V,method,option,
                    % relaxFactor[,dis_opt][,lp_solver]) so an SDPA layer placed
                    % in a composed network does not hit the 'otherwise' error on
                    % every reach call. Only method is consumed (the trailing
                    % reach options are unused by the bounds path); lp_solver
                    % stays at its default.
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

            % Dispatch based on method. Attention is nonlinear (bilinear QK^T +
            % softmax + A*V), so every star-family method routes to the same SOUND
            % over-approximation (there is no exact star reach to fall back to).
            if any(strcmp(method, {'approx-star','exact-star','abs-dom'})) || contains(method, 'relax-star')
                S = obj.reach_star_approx(Q_set, K_set, V_set, lp_solver);
            elseif strcmp(method, 'approx-zono')
                S = obj.reach_zono_approx(Q_set, K_set, V_set);
            else
                error('Unsupported reachability method: %s', method);
            end
        end

        function S = reach_star_approx(obj, Q_set, K_set, V_set, lp_solver)
            % SOUND over-approximate reachability of softmax(scale*Q*K')*V as a
            % Star, valid for MULTI-TOKEN attention. Q,K,V are matrix-Stars whose
            % column-major reshape is [N x D] (D = QueryDim, N = seq length).
            %
            %   single token (N=1): softmax over one key == 1, so attention == V,
            %     and we return V_set unchanged (exact).
            %   multi token: SoftmaxAttn.singleHeadAttn computes sound score bounds
            %     (Rump interval QK'), the exact correlated row-softmax bound, and a
            %     sign-aware symbolic A*V envelope that keeps V's predicates -- so
            %     the output stays correlated with the value path (not a box-lift).
            if nargin < 5
                lp_solver = 'linprog'; %#ok<NASGU>  (bounds use the cheap estimate path)
            end
            D = obj.QueryDim;
            if isempty(D) || D == 0
                error('ScaledDotProductAttentionLayer:noDim', ...
                    'QueryDim must be set to reach (head/key dimension).');
            end
            N = round(Q_set.dim / D);
            if N * D ~= Q_set.dim
                error('ScaledDotProductAttentionLayer:shape', ...
                    'Q dim %d not a multiple of QueryDim %d', Q_set.dim, D);
            end
            if N <= 1
                S = V_set;                       % single token: attention == V (exact)
                return;
            end
            if V_set.dim ~= N*D
                error('ScaledDotProductAttentionLayer:shape', ...
                    'V dim %d ~= N*D %d (this sound path needs ValueDim==QueryDim)', V_set.dim, N*D);
            end
            S = SoftmaxAttn.singleHeadAttn(Q_set, K_set, V_set, obj.Scale, [N D], 'estimate');
        end

        function Z = reach_zono_approx(obj, Q_set, K_set, V_set)
            % Over-approximate reachability returning a (box) Zonotope. Uses the
            % same SOUND multi-token star reach, then encloses it in an axis-aligned
            % zonotope. Q,K,V may be Zono/Star (toStar handles either).
            if isa(Q_set,'Zono'), Q_set = Q_set.toStar(); end
            if isa(K_set,'Zono'), K_set = K_set.toStar(); end
            if isa(V_set,'Zono'), V_set = V_set.toStar(); end
            S = obj.reach_star_approx(Q_set, K_set, V_set, 'linprog');
            [out_lb, out_ub] = S.estimateRanges();
            center = (out_lb + out_ub) / 2;
            generators = diag((out_ub - out_lb) / 2);
            Z = Zono(center, generators);
        end

        function [out_lb, out_ub] = compute_attention_bounds(obj, ~, ~, ~, ~, V_lb, V_ub)
            % Compute bounds on attention output for the SINGLE-TOKEN case.
            %
            % The reach_*_approx routines flatten Q, K, V into vector Stars
            % (treating each as a single-token query/key/value vector). For
            % seq_len = 1 the softmax is taken over a single key, so
            % attn_weights == 1 identically, and Attention(Q,K,V) = V, so
            % returning V's bounds is exact (sound and tight). Q and K unused.
            %
            % SOUNDNESS GUARD (here, not in reach_star_approx, so BOTH the star
            % and zono reach paths are protected -- the zono path used to bypass
            % it, finding [0]): for seq_len > 1 the output is a cross-token convex
            % combination that can leave any single token's value range, so
            % returning V's bounds would be UNSOUND. If the layer carries a
            % per-token ValueDim and the input is an exact multiple of it (>1
            % tokens), refuse. When ValueDim is 0/unset the layer operates under
            % its documented single-token contract (the whole flattened input is
            % ONE token, V-passthrough exact); set ValueDim to arm this guard for
            % multi-token misuse. A sound multi-token bound is future work.
            n_v = numel(V_lb);
            if ~isempty(obj.ValueDim) && obj.ValueDim > 0 && ...
                    n_v > obj.ValueDim && mod(n_v, obj.ValueDim) == 0
                error('ScaledDotProductAttentionLayer:multiTokenUnsound', ...
                    ['Multi-token attention reach (seq_len=%d) is not implemented; the ' ...
                     'single-token bounds path would be UNSOUND. Use a per-token model, ' ...
                     'or implement the per-coordinate cross-token softmax bound ' ...
                     '(see Softmax.compute_softmax_bounds).'], n_v / obj.ValueDim);
            end
            out_lb = V_lb;
            out_ub = V_ub;
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
