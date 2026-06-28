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
        % Operand shapes [m k] and [k n] (out = [m n]). Required to enable a
        % SOUND interval matmul (the operand shapes are not recoverable from the
        % flattened set dimensions). Set them from the importer manifest or when
        % constructing the layer by hand; left empty the layer refuses to reach.
        LeftShape = [];
        RightShape = [];
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

        function S = reach_single_input(obj, inputs)
            % SOUND set@set matrix product via the Rump interval matmul
            % (SoftmaxAttn.bilinearMatMulStar): concretize each operand to its
            % sound per-dimension range, enclose A*B, re-encode as a box Star.
            % Requires LeftShape/RightShape; refuses otherwise (sound-by-refusal:
            % the prior element-wise interval product was NOT a sound matmul bound).
            if isempty(obj.LeftShape) || isempty(obj.RightShape)
                error('DynamicMatmulLayer:reachNeedsShapes', ...
                    ['no operand shapes set: a SOUND interval matmul needs the ' ...
                     '[m k] x [k n] shapes. Set obj.LeftShape/obj.RightShape ' ...
                     '(from the importer manifest), or use a sampling-based method.']);
            end
            A = inputs{1}; B = inputs{2};
            S = SoftmaxAttn.bilinearMatMulStar(A, B, obj.LeftShape, obj.RightShape, 1.0, 'estimate');
        end

        function IS = reach(varargin)
            obj = varargin{1};
            in_images = varargin{2};
            method = varargin{3};
            isStar = any(strcmp(method, {'approx-star','exact-star','abs-dom'})) || contains(method, "relax-star");
            isZono = strcmp(method, 'approx-zono');
            if ~(isStar || isZono)
                error('Unknown reach method for DynamicMatmulLayer: %s', method);
            end
            if strcmp(method, 'exact-star')
                % A product of two dynamic operands is bilinear, so there is no
                % complete exact-star reach; warn so callers do not assume
                % 'exact-star' is tighter than 'approx-star' here (both return the
                % same sound box over-approximation).
                warning('DynamicMatmulLayer:exactStarApprox', ...
                    ['dynamic MatMul is bilinear: no exact star reach exists. ' ...
                     '''exact-star'' returns the SAME sound box over-approximation ' ...
                     'as ''approx-star'' (not complete or tighter).']);
            end
            S = obj.reach_single_input(in_images);     % sound box Star
            if isZono
                [lb, ub] = S.estimateRanges();         % box star -> exact box zono
                IS = Zono((lb+ub)/2, diag((ub-lb)/2));
            else
                IS = S;
            end
        end
    end

    methods
        function obj = toGPU(obj), end
        function obj = changeParamsPrecision(obj, ~), end
    end
end
