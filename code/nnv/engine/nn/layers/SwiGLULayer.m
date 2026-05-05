classdef SwiGLULayer < handle
    % SwiGLULayer - Gated Linear Unit with SiLU activation
    %
    % SwiGLU(x) = SiLU(x * W_gate + b_gate) .* (x * W_up + b_up)
    %
    % This is the gated activation used in modern LLMs:
    %   - SmolLM2
    %   - OLMo
    %   - Llama family
    %   - PaLM
    %
    % The layer computes:
    %   gate = SiLU(x @ W_gate + b_gate)
    %   up = x @ W_up + b_up
    %   output = gate .* up
    %
    % Reference: https://arxiv.org/pdf/2002.05202 (GLU Variants)
    %
    % Author: NNV Team
    % Date: November 2025

    properties
        Name = 'swiglu';
        NumInputs = 1;
        InputNames = {'in'};
        NumOutputs = 1;
        OutputNames = {'out'};

        % Weight matrices
        W_gate = [];    % Gate projection weights [input_dim x hidden_dim]
        b_gate = [];    % Gate projection bias [hidden_dim x 1]
        W_up = [];      % Up projection weights [input_dim x hidden_dim]
        b_up = [];      % Up projection bias [hidden_dim x 1]

        % Dimensions
        InputDim = 0;
        HiddenDim = 0;
    end

    methods

        function obj = SwiGLULayer(varargin)
            % Constructor
            % Usage:
            %   SwiGLULayer(W_gate, W_up)
            %   SwiGLULayer(W_gate, b_gate, W_up, b_up)
            %   SwiGLULayer(W_gate, b_gate, W_up, b_up, name)

            switch nargin
                case 0
                    % Default empty layer
                case 2
                    obj.W_gate = varargin{1};
                    obj.W_up = varargin{2};
                    obj.InputDim = size(obj.W_gate, 1);
                    obj.HiddenDim = size(obj.W_gate, 2);
                    obj.b_gate = zeros(obj.HiddenDim, 1);
                    obj.b_up = zeros(obj.HiddenDim, 1);
                case 4
                    obj.W_gate = varargin{1};
                    obj.b_gate = varargin{2};
                    obj.W_up = varargin{3};
                    obj.b_up = varargin{4};
                    obj.InputDim = size(obj.W_gate, 1);
                    obj.HiddenDim = size(obj.W_gate, 2);
                case 5
                    obj.W_gate = varargin{1};
                    obj.b_gate = varargin{2};
                    obj.W_up = varargin{3};
                    obj.b_up = varargin{4};
                    obj.Name = varargin{5};
                    obj.InputDim = size(obj.W_gate, 1);
                    obj.HiddenDim = size(obj.W_gate, 2);
                otherwise
                    error('Invalid number of arguments');
            end

            % Ensure biases are column vectors
            if ~isempty(obj.b_gate)
                obj.b_gate = obj.b_gate(:);
            end
            if ~isempty(obj.b_up)
                obj.b_up = obj.b_up(:);
            end
        end

    end

    methods  % Evaluation

        function y = evaluate(obj, x)
            % Evaluate SwiGLU: gate .* up
            % @x: input vector [input_dim x 1] or [input_dim x batch]
            % @y: output [hidden_dim x 1] or [hidden_dim x batch]

            x = x(:);  % Ensure column vector

            % Gate path: SiLU(x @ W_gate + b_gate)
            gate_linear = obj.W_gate' * x + obj.b_gate;
            gate = SiLU.evaluate(gate_linear);

            % Up path: x @ W_up + b_up
            up = obj.W_up' * x + obj.b_up;

            % Element-wise multiplication
            y = gate .* up;
        end

    end

    methods  % Reachability

        function R = reach(obj, varargin)
            % Main reachability method
            % @I: input Star or ImageStar
            % @method: 'approx-star', etc.

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
            % Reachability for Star input
            %
            % SwiGLU involves:
            % 1. Two linear projections (gate and up)
            % 2. SiLU on gate path
            % 3. Element-wise multiplication
            %
            % We use a conservative overapproximation based on bounds.

            if ~isa(I, 'Star')
                error('Input must be a Star');
            end

            % Get input bounds
            n = I.dim;
            [in_lb, in_ub] = I.getRanges;

            % Step 1: Compute bounds on gate_linear = W_gate' * x + b_gate
            gate_linear_lb = zeros(obj.HiddenDim, 1);
            gate_linear_ub = zeros(obj.HiddenDim, 1);

            for i = 1:obj.HiddenDim
                w = obj.W_gate(:, i);
                b = obj.b_gate(i);

                % Interval arithmetic for linear combination
                pos_w = max(w, 0);
                neg_w = min(w, 0);

                gate_linear_lb(i) = pos_w' * in_lb + neg_w' * in_ub + b;
                gate_linear_ub(i) = pos_w' * in_ub + neg_w' * in_lb + b;
            end

            % Step 2: Compute bounds on gate = SiLU(gate_linear)
            gate_lb = zeros(obj.HiddenDim, 1);
            gate_ub = zeros(obj.HiddenDim, 1);

            for i = 1:obj.HiddenDim
                l = gate_linear_lb(i);
                u = gate_linear_ub(i);

                % SiLU bounds on interval [l, u]
                [gate_lb(i), gate_ub(i)] = SwiGLULayer.silu_bounds(l, u);
            end

            % Step 3: Compute bounds on up = W_up' * x + b_up
            up_lb = zeros(obj.HiddenDim, 1);
            up_ub = zeros(obj.HiddenDim, 1);

            for i = 1:obj.HiddenDim
                w = obj.W_up(:, i);
                b = obj.b_up(i);

                pos_w = max(w, 0);
                neg_w = min(w, 0);

                up_lb(i) = pos_w' * in_lb + neg_w' * in_ub + b;
                up_ub(i) = pos_w' * in_ub + neg_w' * in_lb + b;
            end

            % Step 4: Compute bounds on output = gate .* up
            % Using interval multiplication
            out_lb = zeros(obj.HiddenDim, 1);
            out_ub = zeros(obj.HiddenDim, 1);

            for i = 1:obj.HiddenDim
                [out_lb(i), out_ub(i)] = SwiGLULayer.interval_multiply(...
                    gate_lb(i), gate_ub(i), up_lb(i), up_ub(i));
            end

            % Create output Star from bounds
            R = Star(out_lb, out_ub);
        end

        function R = reach_imagestar(obj, I, method)
            % Reachability for ImageStar input

            if ~isa(I, 'ImageStar')
                error('Input must be an ImageStar');
            end

            % Convert to Star, process, return as Star
            % (SwiGLU typically receives flattened input)
            S = I.toStar;
            R = obj.reach_star(S, method);
        end

        function R = reach_star_single_input(obj, I, method, ~, ~, ~)
            % Compatibility method for NNV layer interface
            R = obj.reach_star(I, method);
        end

    end

    methods(Static)

        function [lb, ub] = silu_bounds(l, u)
            % Compute bounds on SiLU(x) for x in [l, u]

            if l == u
                lb = SiLU.evaluate(l);
                ub = lb;
                return;
            end

            % Evaluate at endpoints
            y_l = SiLU.evaluate(l);
            y_u = SiLU.evaluate(u);

            % SiLU minimum is at x ≈ -1.278 with value ≈ -0.278
            x_min = -1.2784645427610737;
            y_min = -0.27846454276107373;

            % Check if minimum is in interval
            if l <= x_min && x_min <= u
                lb = y_min;
            else
                lb = min(y_l, y_u);
            end

            ub = max(y_l, y_u);
        end

        function [lb, ub] = interval_multiply(a_lb, a_ub, b_lb, b_ub)
            % Compute bounds on a * b where a in [a_lb, a_ub], b in [b_lb, b_ub]

            products = [a_lb * b_lb, a_lb * b_ub, a_ub * b_lb, a_ub * b_ub];
            lb = min(products);
            ub = max(products);
        end

        function L = parse(layer)
            % Parse a MATLAB/ONNX SwiGLU layer
            % Note: SwiGLU is typically not a single layer in frameworks
            % This is a placeholder for custom parsing

            error('SwiGLU parsing not implemented - use constructor directly');
        end

    end

    methods  % Helper methods

        function obj = set_weights(obj, W_gate, b_gate, W_up, b_up)
            % Set layer weights
            obj.W_gate = W_gate;
            obj.b_gate = b_gate(:);
            obj.W_up = W_up;
            obj.b_up = b_up(:);
            obj.InputDim = size(W_gate, 1);
            obj.HiddenDim = size(W_gate, 2);
        end

        function obj = toGPU(obj)
            % Move weights to GPU
            obj.W_gate = gpuArray(obj.W_gate);
            obj.b_gate = gpuArray(obj.b_gate);
            obj.W_up = gpuArray(obj.W_up);
            obj.b_up = gpuArray(obj.b_up);
        end

        function obj = changeParamsPrecision(obj, precision)
            % Change parameter precision
            if strcmp(precision, 'single')
                obj.W_gate = single(obj.W_gate);
                obj.b_gate = single(obj.b_gate);
                obj.W_up = single(obj.W_up);
                obj.b_up = single(obj.b_up);
            elseif strcmp(precision, 'double')
                obj.W_gate = double(obj.W_gate);
                obj.b_gate = double(obj.b_gate);
                obj.W_up = double(obj.W_up);
                obj.b_up = double(obj.b_up);
            end
        end

    end

end
