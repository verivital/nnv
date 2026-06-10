classdef DynamicMatmulLayer < handle
    % Matrix product of two dynamic tensors: out = inputs{1} * inputs{2}.
    % Used for ONNX MatMul ops where neither operand is an initializer
    % (e.g. attention Q*K^T and attn_weights*V in transformers).
    %
    % evaluate is exact (matrix multiplication).
    % reach is NOT implemented soundly and therefore ERRORS: the layer has no
    % operand-shape metadata, and the previous element-wise interval
    % approximation was NOT a sound bound for a matrix product
    % (out(i,j) = sum_k A(i,k)*B(k,j) sums k bilinear terms and changes
    % shape). A sound interval matmul needs the [m x k], [k x n] shapes --
    % wire them through the importer manifest before enabling reach (TODO).

    properties
        Name = 'DynamicMatmulLayer';
        NumInputs = 2;
        InputNames = {'in1', 'in2'};
        NumOutputs = 1;
        OutputNames = {'out'};
    end

    methods
        function obj = DynamicMatmulLayer(varargin)
            switch nargin
                case 5
                    obj.Name        = varargin{1};
                    obj.NumInputs   = varargin{2};
                    obj.NumOutputs  = varargin{3};
                    obj.InputNames  = varargin{4};
                    obj.OutputNames = varargin{5};
                case 1
                    obj.Name = varargin{1};
                otherwise
                    error('DynamicMatmulLayer: ctor takes 1 or 5 args');
            end
            if obj.NumInputs ~= 2
                error('DynamicMatmulLayer requires exactly 2 inputs');
            end
        end

        function out = evaluate(~, inputs)
            % Standard matrix multiplication. ONNX MatMul broadcasts on
            % leading dimensions, so for rank-3+ tensors we apply matmul
            % along the last two dims per leading slice.
            A = inputs{1};
            B = inputs{2};
            if ndims(A) <= 2 && ndims(B) <= 2
                out = A * B;
                return;
            end
            % Broadcasting matmul along last two axes
            sa = size(A); sb = size(B);
            k_a = sa(end-1); k_a_in = sa(end);
            k_b = sb(end-1); k_b_in = sb(end);
            if k_a_in ~= k_b
                error('DynamicMatmulLayer: inner dims mismatch (%d vs %d)', k_a_in, k_b);
            end
            % Reshape to [batch, k_a, k_a_in] and [batch, k_b, k_b_in]
            % then matmul slice-by-slice.
            lead_a = sa(1:end-2);
            lead_b = sb(1:end-2);
            % [6] The slice-by-slice loop below pairs A's i-th leading slice with
            % B's i-th -- correct ONLY when the leading/batch dims are identical.
            % The old code assumed that without checking (`lead = lead_a`), so for
            % e.g. A [2,3,m,k] and B [3,2,k,n] (prod leads equal but ORDER
            % different) it silently multiplied mismatched slices, and for a
            % genuine ONNX leading-dim broadcast (e.g. B [3,k,n] vs A [2,3,m,k])
            % the B reshape just crashed. Require equal leading dims; refuse
            % otherwise rather than return a silently-wrong product. (Full
            % numpy-style leading-dim broadcasting is future work.)
            if ~isequal(lead_a, lead_b)
                error('DynamicMatmulLayer:leadingDimsDiffer', ...
                    ['leading/batch dims differ (%s vs %s); only equal leading ' ...
                     'dims are supported. ONNX-style leading-dim broadcasting is ' ...
                     'not implemented -- refusing to pair mismatched slices.'], ...
                    mat2str(lead_a), mat2str(lead_b));
            end
            lead = lead_a;
            n_lead = max(prod(lead), 1);
            A_flat = reshape(A, [n_lead, k_a, k_a_in]);
            B_flat = reshape(B, [n_lead, k_b, k_b_in]);
            out_flat = zeros(n_lead, k_a, k_b_in, 'like', A);
            for i = 1:n_lead
                Ai = squeeze(A_flat(i, :, :));
                Bi = squeeze(B_flat(i, :, :));
                if isvector(Ai), Ai = reshape(Ai, k_a, k_a_in); end
                if isvector(Bi), Bi = reshape(Bi, k_b, k_b_in); end
                out_flat(i, :, :) = Ai * Bi;
            end
            out = reshape(out_flat, [lead, k_a, k_b_in]);
        end

        function out = reach_single_input(~, inputs) %#ok<INUSD,STOUT>
            % SOUNDNESS: refuse rather than approximate wrongly. The previous
            % implementation took element-wise interval products of the two
            % flattened operands -- that is NOT a sound bound for a matrix
            % product (each out(i,j) sums k bilinear terms and the output
            % shape differs from either input), and its fallbacks silently
            % passed an operand through as "the reachable set". Without the
            % operand shapes a sound interval matmul cannot be computed here.
            % evaluate() remains exact.
            error('DynamicMatmulLayer:reachNotImplemented', ...
                ['no SOUND reachability is implemented for dynamic MatMul ' ...
                 '(needs operand shapes for a per-element sum-of-bilinear-terms ' ...
                 'interval bound). Refusing to return an unsound set; use a ' ...
                 'sampling-based method (e.g. cp-star) for networks with ' ...
                 'dynamic MatMul, or wire shapes through the importer.']);
        end

        function IS = reach(varargin)
            obj = varargin{1};
            in_images = varargin{2};
            method = varargin{3};
            if any(strcmp(method, {'approx-star','exact-star','abs-dom','approx-zono'})) || contains(method, "relax-star")
                IS = obj.reach_single_input(in_images);
            else
                error('Unknown reach method for DynamicMatmulLayer: %s', method);
            end
        end
    end

    methods
        function obj = toGPU(obj), end
        function obj = changeParamsPrecision(obj, ~), end
    end
end
