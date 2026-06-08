classdef ElementwiseDivisionLayer < handle
    % Element-wise division layer (a ./ b).
    % Used for ONNX Div ops on dynamic tensors.
    %
    % Soundness note: reach uses an interval over-approximation. Division
    % by an interval that straddles zero is not finite — those cases are
    % flagged with a warning and the output is passed through.

    properties
        Name = 'DivLayer';
        NumInputs = 2;
        InputNames = {'in1', 'in2'};
        NumOutputs = 1;
        OutputNames = {'out'};
    end

    methods
        function obj = ElementwiseDivisionLayer(varargin)
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
                    error('ElementwiseDivisionLayer: ctor takes 1 or 5 args');
            end
            if obj.NumInputs ~= 2
                error('ElementwiseDivisionLayer requires exactly 2 inputs');
            end
        end

        function out = evaluate(~, inputs)
            out = inputs{1} ./ inputs{2};
        end

        function out = reach_single_input(~, inputs)
            if ~iscell(inputs) || length(inputs) ~= 2
                out = inputs;
                return;
            end
            try
                [a_lb, a_ub] = inputs{1}.getRanges();
                [b_lb, b_ub] = inputs{2}.getRanges();
                if any(b_lb <= 0 & b_ub >= 0)
                    warning('ElementwiseDivisionLayer: divisor straddles zero, reach unsound');
                    out = inputs{1};
                    return;
                end
                a_lb = a_lb(:); a_ub = a_ub(:); b_lb = b_lb(:); b_ub = b_ub(:);
                c = [a_lb./b_lb, a_lb./b_ub, a_ub./b_lb, a_ub./b_ub];
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
                error('Unknown reach method for ElementwiseDivisionLayer: %s', method);
            end
        end
    end

    methods
        function obj = toGPU(obj), end
        function obj = changeParamsPrecision(obj, ~), end
    end
end
