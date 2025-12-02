classdef SiLULayer < ActivationFunctionLayer
    % SiLULayer - Sigmoid Linear Unit (SiLU/Swish) activation layer
    %
    % SiLU(x) = x * sigmoid(x) = x / (1 + exp(-x))
    %
    % This activation is used in modern LLMs including:
    %   - SmolLM2
    %   - OLMo
    %   - Llama family
    %   - GPT-NeoX
    %
    % Properties:
    %   - Smooth and differentiable everywhere
    %   - Non-monotonic (has a minimum at x ≈ -1.278)
    %   - Bounded below: min(SiLU(x)) ≈ -0.278
    %   - Self-gated: output depends on both x and sigmoid(x)
    %
    % Reference: https://arxiv.org/abs/1710.05941
    %
    % Author: NNV Team
    % Date: November 2025

    methods

        function obj = SiLULayer(varargin)
            % Constructor
            % Usage:
            %   SiLULayer()
            %   SiLULayer(name)
            %   SiLULayer(name, numInputs, inputNames, numOutputs, outputNames)

            obj = obj@ActivationFunctionLayer(varargin);

            % Set default name if not provided
            if nargin == 0
                obj.Name = 'silu';
            end
        end

    end

    methods  % Evaluation method

        function y = evaluate(~, input)
            % Evaluate SiLU activation
            % @input: N-dimensional array
            % @y: output with same shape as input

            n = size(input);
            N = prod(n);

            I = reshape(input, [N, 1]);
            y = SiLU.evaluate(I);
            y = reshape(y, n);
        end

    end

    methods  % Reachability methods

        function images = reach_star_single_input(~, in_image, method, relaxFactor, dis_opt, lp_solver)
            % Reachability using ImageStar or Star
            % @in_image: ImageStar or Star input set
            % @method: 'exact-star', 'approx-star', or 'abs-dom'
            % @relaxFactor: for approx-star method
            % @images: ImageStar or Star output set(s)

            if ~isa(in_image, 'ImageStar') && ~isa(in_image, 'Star')
                error('Input must be an ImageStar or Star');
            end

            if isa(in_image, 'ImageStar')
                h = in_image.height;
                w = in_image.width;
                c = in_image.numChannel;

                % Convert to Star, compute reachability, convert back
                Y = SiLU.reach(in_image.toStar, method, [], relaxFactor, dis_opt, lp_solver);

                n = length(Y);
                images(n) = ImageStar;
                for i = 1:n
                    images(i) = Y(i).toImageStar(h, w, c);
                end
            else
                images = SiLU.reach(in_image, method, [], relaxFactor, dis_opt, lp_solver);
            end
        end

        function image = reach_zono(~, in_image)
            % Reachability using ImageZono or Zono
            % @in_image: ImageZono or Zono input set

            if ~isa(in_image, 'ImageZono') && ~isa(in_image, 'Zono')
                error('Input must be an ImageZono or Zono');
            end

            if isa(in_image, 'ImageZono')
                h = in_image.height;
                w = in_image.width;
                c = in_image.numChannels;
                In = in_image.toZono;
                Y = SiLU.reach(In, 'approx-zono');
                image = Y.toImageZono(h, w, c);
            else
                image = SiLU.reach(in_image, 'approx-zono');
            end
        end

    end

    methods(Static)

        function L = parse(layer)
            % Parse a MATLAB/ONNX SiLU/Swish layer
            % @layer: MATLAB or ONNX layer object
            % @L: SiLULayer object

            % Check for known SiLU/Swish layer types
            valid_types = {...
                'nnet.onnx.layer.SwishLayer', ...
                'nnet.keras.layer.SwishLayer', ...
                'SwishLayer', ...
                'nnet.cnn.layer.SwishLayer' ...
            };

            is_valid = false;
            for i = 1:length(valid_types)
                if isa(layer, valid_types{i})
                    is_valid = true;
                    break;
                end
            end

            % Also check by name/type string
            if ~is_valid
                layer_class = class(layer);
                if contains(lower(layer_class), 'swish') || contains(lower(layer_class), 'silu')
                    is_valid = true;
                end
            end

            if ~is_valid
                error('Input is not a recognized SiLU/Swish layer type');
            end

            % Create layer with properties from parsed layer
            if isprop(layer, 'Name')
                name = layer.Name;
            else
                name = 'silu';
            end

            if isprop(layer, 'NumInputs')
                numInputs = layer.NumInputs;
            else
                numInputs = 1;
            end

            if isprop(layer, 'InputNames')
                inputNames = layer.InputNames;
            else
                inputNames = {'in'};
            end

            if isprop(layer, 'NumOutputs')
                numOutputs = layer.NumOutputs;
            else
                numOutputs = 1;
            end

            if isprop(layer, 'OutputNames')
                outputNames = layer.OutputNames;
            else
                outputNames = {'out'};
            end

            L = SiLULayer(name, numInputs, inputNames, numOutputs, outputNames);
        end

    end

end
