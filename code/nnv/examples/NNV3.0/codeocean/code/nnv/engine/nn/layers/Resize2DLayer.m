classdef Resize2DLayer < handle
    % Resize2DLayer for NNV
    % Implements resizing ImageStars by OutputSize or Scale
    
    properties
        Name = 'Resize2DLayer';
        NumInputs = 1;
        InputNames = {'in1'};
        NumOutputs = 1;
        OutputNames = {'out'};
        
        OutputSize;            % [H, W] or []
        Scale;                 % [scaleH, scaleW] or []
        Method;                % 'nearest', 'bilinear', etc.
        Antialiasing = false;  % true/false
        GeometricTransformMode = 'half-pixel'; % placeholder
        NearestRoundingMode = 'round';         % placeholder
    end
    
    methods
        function obj = Resize2DLayer(varargin)
            % Constructor: expect (Name, NumInputs, NumOutputs, InputNames, OutputNames, OutputSize, Scale, Method, Antialiasing, GeometricTransformMode, NearestRoundingMode)
            if nargin == 11
                obj.Name = varargin{1};
                obj.NumInputs = varargin{2};
                obj.NumOutputs = varargin{3};
                obj.InputNames = varargin{4};
                obj.OutputNames = varargin{5};
                obj.OutputSize = varargin{6};
                obj.Scale = varargin{7};
                obj.Method = varargin{8};
                obj.Antialiasing = varargin{9};
                obj.GeometricTransformMode = varargin{10};
                obj.NearestRoundingMode = varargin{11};
            else
                error('Invalid number of input arguments');
            end
        end
        
        % Reachability with a single ImageStar
        function IS = reach_single_input(obj, input)
            if ~isa(input, 'ImageStar')
                error('Input must be an ImageStar');
            end

            % Determine target size
            if ~isempty(obj.OutputSize)
                targetSize = obj.OutputSize;
            elseif ~isempty(obj.Scale)
                inputSize = size(input.V(:,:,1,1)); % [H, W]
                targetSize = round(obj.Scale .* inputSize);
            else
                error('Resize2DLayer must specify either OutputSize or Scale.');
            end

            numPred = size(input.V, 4);
            outV = cell(1, numPred);

            for i = 1:numPred
                outV{i} = imresize(input.V(:,:,:,i), targetSize, ...
                    'Method', obj.Method, ...
                    'Antialiasing', obj.Antialiasing);
            end

            Vout = cat(4, outV{:});
            IS = ImageStar(Vout, input.C, input.d, input.pred_lb, input.pred_ub);
        end
        
        % Reachability with multiple inputs
        function S = reach_multipleInputs(obj, inputs, option)
            n = length(inputs);
            if isa(inputs(1), 'ImageStar')
                S(n) = ImageStar;
            else
                error('Unknown input data set');
            end

            if strcmp(option, 'parallel')
                parfor i = 1:n
                    S(i) = obj.reach_single_input(inputs(i));
                end
            elseif strcmp(option, 'single') || isempty(option)
                for i = 1:n
                    S(i) = obj.reach_single_input(inputs(i));
                end
            else
                error('Unknown computation option');
            end
        end
        
        % Main reach method dispatcher
        function IS = reach(obj, varargin)
            narginchk(2,4)
            in_image = varargin{1};
            if length(varargin) >= 3
                option = varargin{3};
            else
                option = [];
            end

            if strcmp(option, 'parallel') || strcmp(option, 'single') || isempty(option)
                if iscell(in_image) || length(in_image) > 1
                    IS = obj.reach_multipleInputs(in_image, option);
                else
                    IS = obj.reach_single_input(in_image);
                end
            else
                error('Unknown computation option');
            end
        end
    end
    
    methods (Static)
        function L = parse(layer)
            if ~isa(layer, 'nnet.cnn.layer.Resize2DLayer')
                error('Input is not a Resize2DLayer');
            end

            % Extract everything
            name = layer.Name;
            numInputs = layer.NumInputs;
            numOutputs = layer.NumOutputs;
            inputNames = layer.InputNames;
            outputNames = layer.OutputNames;
            outputSize = [];
            scale = [];
            if isprop(layer, 'OutputSize')
                outputSize = layer.OutputSize;
            end
            if isprop(layer, 'Scale')
                scale = layer.Scale;
            end
            method = layer.Method;
            antialiasing = false;
            if isprop(layer, 'Antialiasing')
                antialiasing = layer.Antialiasing;
            end
            geomMode = 'half-pixel';
            if isprop(layer, 'GeometricTransformMode')
                geomMode = layer.GeometricTransformMode;
            end
            nearestMode = 'round';
            if isprop(layer, 'NearestRoundingMode')
                nearestMode = layer.NearestRoundingMode;
            end

            L = Resize2DLayer(name, numInputs, numOutputs, inputNames, outputNames, ...
                              outputSize, scale, method, antialiasing, geomMode, nearestMode);
        end
    end
end



