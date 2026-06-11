classdef ElementwiseProductLayer < handle
    % Element-wise product (Hadamard) layer.
    % Multiplies two or more dynamic input tensors element-wise.
    %
    % Used for ONNX Mul ops where both operands are dynamic (intermediate
    % tensors), e.g. Lyapunov function bilinear forms `x .* (M*x)` and
    % attention scaling `scores .* scale_tensor`.
    %
    % Soundness note: reach_single_input uses an interval over-approximation
    % of the bilinear product. Tighter Star-set bounds would require a custom
    % bilinear abstract transformer (TODO).

    properties
        Name = 'ProductLayer';
        NumInputs = 2;
        InputNames = {'in1', 'in2'};
        NumOutputs = 1;
        OutputNames = {'out'};
    end

    methods
        function obj = ElementwiseProductLayer(varargin)
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
                    error('ElementwiseProductLayer: ctor takes 1 or 5 args');
            end
            if obj.NumInputs < 2
                error('ElementwiseProductLayer needs >= 2 inputs');
            end
        end

        function out = evaluate(~, inputs)
            % @inputs: cell array of arrays of equal/broadcastable shape.
            % @out: their element-wise product.
            out = inputs{1};
            for k = 2:length(inputs)
                out = out .* inputs{k};
            end
        end

        function out = reach_single_input(~, inputs)
            % Interval over-approximation of the elementwise product (sound:
            % the four corner products bound a bilinear term over a box).
            % For a constant scale use ElementwiseAffineLayer (linear, exact).
            %
            % SOUNDNESS: never silently pass an operand through as "the
            % reachable set" -- all failure modes raise instead (the previous
            % catch->passthrough was explicitly unsound).
            if ~iscell(inputs) || length(inputs) < 2
                error('ElementwiseProductLayer:badInputs', ...
                    'reach requires a cell of >= 2 input sets; got %s', class(inputs));
            end
            lb = []; ub = [];
            for k = 1:length(inputs)
                [lk, uk] = inputs{k}.getRanges();
                lk = lk(:); uk = uk(:);
                if isempty(lb)
                    lb = lk; ub = uk;
                else
                    % element-wise interval product
                    c = [lb.*lk, lb.*uk, ub.*lk, ub.*uk];
                    lb = min(c, [], 2);
                    ub = max(c, [], 2);
                end
            end
            out = Star(lb, ub);
        end

        function S = reach_multipleInputs(obj, inputs, ~)
            % ElementwiseProduct multiplies its operands TOGETHER (it is a
            % binary/N-ary op, NOT a per-element batch), so "multiple inputs"
            % means the cell of operands -- exactly what reach_single_input
            % consumes. [47]: the previous body forwarded inputs(i) (a single
            % non-cell set) into reach_single_input, whose >=2-element-cell guard
            % made EVERY call error -- a dead path that also dropped all operands
            % but the indexed one. Delegate to the operand handler instead.
            S = obj.reach_single_input(inputs);
        end

        function IS = reach(varargin)
            obj = varargin{1};
            in_images = varargin{2};
            method = varargin{3};
            if any(strcmp(method, {'approx-star','exact-star','abs-dom','approx-zono'})) || contains(method, "relax-star")
                IS = obj.reach_single_input(in_images);
            else
                error('Unknown reach method for ElementwiseProductLayer: %s', method);
            end
        end
    end

    methods
        function obj = toGPU(obj), end
        function obj = changeParamsPrecision(obj, ~), end
    end
end
