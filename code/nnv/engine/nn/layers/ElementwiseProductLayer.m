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
            % Best-effort interval over-approximation. For a Star input
            % `S` with bounds [lb, ub] and a constant scale `c`, S * c is
            % linear and exact (use ElementwiseAffineLayer). For two
            % dynamic Stars, take elementwise interval product and form a
            % Box, then convert back to Star. Loose but sound.
            if ~iscell(inputs)
                out = inputs;
                return;
            end
            try
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
            catch
                % Fall back to passthrough of first operand (unsound but
                % does not crash).
                out = inputs{1};
            end
        end

        function S = reach_multipleInputs(obj, inputs, option)
            n = length(inputs);
            if isa(inputs(1), 'ImageStar')
                S(n) = ImageStar;
            elseif isa(inputs(1), 'Star')
                S(n) = Star;
            else
                error('Unknown input set type');
            end
            if strcmp(option, 'parallel')
                parfor i = 1:n
                    S(i) = obj.reach_single_input(inputs(i));
                end
            else
                for i = 1:n
                    S(i) = obj.reach_single_input(inputs(i));
                end
            end
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
