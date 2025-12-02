classdef RMSNormLayer < handle
    % RMSNormLayer - Root Mean Square Layer Normalization
    %
    % RMSNorm(x) = x / sqrt(mean(x^2) + eps) * gamma
    %
    % Unlike LayerNorm, RMSNorm does NOT subtract the mean (no centering),
    % only scales by the root mean square. This is used in modern LLMs:
    %   - SmolLM2
    %   - OLMo
    %   - Llama family
    %   - GPT-NeoX
    %
    % Reference: https://arxiv.org/abs/1910.07467
    %
    % Author: NNV Team
    % Date: November 2025

    properties
        Name = 'rmsnorm';
        NumInputs = 1;
        InputNames = {'in'};
        NumOutputs = 1;
        OutputNames = {'out'};

        % Parameters
        Gamma = [];     % Scale parameter (learned) [dim x 1]
        Epsilon = 1e-6; % Small constant for numerical stability

        % Dimension
        Dim = 0;
    end

    methods

        function obj = RMSNormLayer(varargin)
            % Constructor
            % Usage:
            %   RMSNormLayer(gamma)
            %   RMSNormLayer(gamma, epsilon)
            %   RMSNormLayer(gamma, epsilon, name)
            %   RMSNormLayer(dim) - creates with ones for gamma

            switch nargin
                case 0
                    % Default
                case 1
                    if isscalar(varargin{1})
                        % Dimension provided
                        obj.Dim = varargin{1};
                        obj.Gamma = ones(obj.Dim, 1);
                    else
                        % Gamma provided
                        obj.Gamma = varargin{1}(:);
                        obj.Dim = length(obj.Gamma);
                    end
                case 2
                    obj.Gamma = varargin{1}(:);
                    obj.Epsilon = varargin{2};
                    obj.Dim = length(obj.Gamma);
                case 3
                    obj.Gamma = varargin{1}(:);
                    obj.Epsilon = varargin{2};
                    obj.Name = varargin{3};
                    obj.Dim = length(obj.Gamma);
                otherwise
                    error('Invalid number of arguments');
            end
        end

    end

    methods  % Evaluation

        function y = evaluate(obj, x)
            % Evaluate RMSNorm
            % @x: input vector [dim x 1] or array
            % @y: normalized output

            x = x(:);  % Ensure column vector

            % Compute RMS = sqrt(mean(x^2))
            rms = sqrt(mean(x.^2) + obj.Epsilon);

            % Normalize and scale
            y = (x / rms) .* obj.Gamma;
        end

    end

    methods  % Reachability

        function R = reach(obj, varargin)
            % Main reachability method

            switch nargin
                case 2
                    I = varargin{1};
                    method = 'approx-star';
                case 3
                    I = varargin{1};
                    method = varargin{2};
                otherwise
                    I = varargin{1};
                    method = varargin{2};
            end

            if isa(I, 'ImageStar')
                R = obj.reach_imagestar(I, method);
            elseif isa(I, 'Star')
                R = obj.reach_star(I, method);
            else
                error('Input must be Star or ImageStar');
            end
        end

        function R = reach_star(obj, I, method)
            % Reachability for Star input using interval bounds
            %
            % RMSNorm is challenging because:
            % 1. x^2 introduces nonlinearity
            % 2. Division by sqrt(sum) is another nonlinearity
            %
            % We use conservative interval arithmetic.

            if ~isa(I, 'Star')
                error('Input must be a Star');
            end

            n = I.dim;

            % Get input bounds
            [in_lb, in_ub] = I.getRanges;

            % Step 1: Bound x^2 element-wise
            x2_lb = zeros(n, 1);
            x2_ub = zeros(n, 1);

            for i = 1:n
                l = in_lb(i);
                u = in_ub(i);

                % x^2 on interval [l, u]
                if l >= 0
                    % Both positive: min at l, max at u
                    x2_lb(i) = l^2;
                    x2_ub(i) = u^2;
                elseif u <= 0
                    % Both negative: min at u, max at l
                    x2_lb(i) = u^2;
                    x2_ub(i) = l^2;
                else
                    % Crosses zero: min is 0, max at further endpoint
                    x2_lb(i) = 0;
                    x2_ub(i) = max(l^2, u^2);
                end
            end

            % Step 2: Bound mean(x^2)
            mean_x2_lb = sum(x2_lb) / n;
            mean_x2_ub = sum(x2_ub) / n;

            % Step 3: Bound sqrt(mean(x^2) + eps)
            rms_lb = sqrt(mean_x2_lb + obj.Epsilon);
            rms_ub = sqrt(mean_x2_ub + obj.Epsilon);

            % Step 4: Bound 1/rms (inverse bounds swap)
            inv_rms_lb = 1 / rms_ub;
            inv_rms_ub = 1 / rms_lb;

            % Step 5: Bound x / rms for each element
            out_lb = zeros(n, 1);
            out_ub = zeros(n, 1);

            for i = 1:n
                % Interval multiplication: x[i] * (1/rms)
                [prod_lb, prod_ub] = RMSNormLayer.interval_multiply(...
                    in_lb(i), in_ub(i), inv_rms_lb, inv_rms_ub);

                % Apply gamma scaling
                if obj.Gamma(i) >= 0
                    out_lb(i) = prod_lb * obj.Gamma(i);
                    out_ub(i) = prod_ub * obj.Gamma(i);
                else
                    out_lb(i) = prod_ub * obj.Gamma(i);
                    out_ub(i) = prod_lb * obj.Gamma(i);
                end
            end

            % Create output Star from bounds
            R = Star(out_lb, out_ub);
        end

        function R = reach_imagestar(obj, I, method)
            % Reachability for ImageStar input

            if ~isa(I, 'ImageStar')
                error('Input must be an ImageStar');
            end

            % Flatten, process, return as Star
            S = I.toStar;
            R = obj.reach_star(S, method);
        end

        function R = reach_star_single_input(obj, I, method, ~, ~, ~)
            % Compatibility method
            R = obj.reach_star(I, method);
        end

    end

    methods(Static)

        function [lb, ub] = interval_multiply(a_lb, a_ub, b_lb, b_ub)
            % Compute bounds on a * b

            products = [a_lb * b_lb, a_lb * b_ub, a_ub * b_lb, a_ub * b_ub];
            lb = min(products);
            ub = max(products);
        end

        function L = parse(layer)
            % Parse a MATLAB/ONNX RMSNorm layer
            % Note: RMSNorm may not be directly available in frameworks

            if isprop(layer, 'Weights') || isprop(layer, 'Scale')
                if isprop(layer, 'Weights')
                    gamma = layer.Weights;
                else
                    gamma = layer.Scale;
                end

                if isprop(layer, 'Epsilon')
                    eps = layer.Epsilon;
                else
                    eps = 1e-6;
                end

                if isprop(layer, 'Name')
                    name = layer.Name;
                else
                    name = 'rmsnorm';
                end

                L = RMSNormLayer(gamma, eps, name);
            else
                error('Cannot parse layer - missing expected properties');
            end
        end

    end

    methods  % Helper methods

        function obj = set_gamma(obj, gamma)
            % Set gamma parameter
            obj.Gamma = gamma(:);
            obj.Dim = length(obj.Gamma);
        end

        function obj = toGPU(obj)
            % Move parameters to GPU
            obj.Gamma = gpuArray(obj.Gamma);
        end

        function obj = changeParamsPrecision(obj, precision)
            % Change parameter precision
            if strcmp(precision, 'single')
                obj.Gamma = single(obj.Gamma);
            elseif strcmp(precision, 'double')
                obj.Gamma = double(obj.Gamma);
            end
        end

    end

end
