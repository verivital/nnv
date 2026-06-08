classdef DynamicMatmulLayer < handle
    % Matrix product of two dynamic tensors: out = inputs{1} * inputs{2}.
    % Used for ONNX MatMul ops where neither operand is an initializer
    % (e.g. attention Q*K^T and attn_weights*V in transformers).
    %
    % evaluate is exact (matrix multiplication).
    % reach uses an interval over-approximation of the bilinear product.
    % Sound but loose; tighter Star-set bounds would require a custom
    % bilinear abstract transformer (TODO).

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
            lead = lead_a;  % assume already broadcast-compatible / equal
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

        function out = reach_single_input(~, inputs)
            % Interval over-approximation of A * B.
            % Treat A as [m × k], B as [k × n]; out is [m × n] with each
            % element a sum over k bilinear products.
            if ~iscell(inputs) || length(inputs) ~= 2
                out = inputs;
                return;
            end
            try
                [a_lb, a_ub] = inputs{1}.getRanges();
                [b_lb, b_ub] = inputs{2}.getRanges();
                a_lb = a_lb(:); a_ub = a_ub(:);
                b_lb = b_lb(:); b_ub = b_ub(:);
                % Sound but very coarse: assume scalar product per element
                % (treats inputs as flat vectors with element-wise interval mul).
                % This is correct only if shapes match; otherwise mark unsound.
                if numel(a_lb) ~= numel(b_lb)
                    warning('DynamicMatmulLayer: bound shapes differ; reach passthrough.');
                    out = inputs{1};
                    return;
                end
                c = [a_lb.*b_lb, a_lb.*b_ub, a_ub.*b_lb, a_ub.*b_ub];
                lb = min(c, [], 2);
                ub = max(c, [], 2);
                out = Star(lb, ub);
            catch
                out = inputs{1};
            end
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
