classdef TransposedConv2DLayer < handle
    % The transpose convolutional 2D layer class
    %   Contain constructor and reachability analysis methods
    % Matlab reference: 
    %   1) https://www.mathworks.com/help/deeplearning/ref/nnet.cnn.layer.transposedconvolution2dlayer.html
    % Deconvolution explanation
    %   1) https://towardsdatascience.com/review-deconvnet-unpooling-layer-semantic-segmentation-55cf8a6e380e
    
    %   Dung Tran: 4/22/2020
    
    properties
        Name = 'transpose_conv_layer';
        % Hyperparameters
        FilterSize = []; % height and width of filters
        NumChannels = 'auto';
        NumFilters = 0; % number of filters
        Stride = [1 1]; % step size for traversing input
        
        CroppingMode = 'manual';
        CroppingSize = [0 0 0 0]; % Outputsize reduction
        
        % Learnable Parmeters/ Used for reachability analysis
        Weights = [];
        Bias = [];
        
        NumInputs = 1;
        InputNames = {'in'};
        NumOutputs = 1;
        OutputNames = {'out'};
    end
    
    
    % setting hyperparameters method
    methods
        
        % constructor of the class
        function obj = TransposedConv2DLayer(varargin)           
            
            switch nargin
                
                case 9 % used for parsing a Matlab conv2dlayer 
                    
                    obj.Name = varargin{1}; 
                    obj.Weights = varargin{2};
                    obj.Bias = varargin{3};
                    obj.CroppingSize = varargin{4};
                    obj.Stride = varargin{5};
                    obj.NumInputs = varargin{6};
                    obj.InputNames = varargin{7};
                    obj.NumOutputs = varargin{8};
                    obj.OutputNames = varargin{9};
                                        
                    w = size(obj.Weights);
                    b = size(obj.Bias);
                    
                    if length(w) ~= 4
                        error('Invalid weights array');
                    end
                    if (w(3) > 1 && length(b) ~= 3) || (w(3) == 1 && length(obj.Bias) ~= 1)
                        error('Invalid biases array');
                    end

                    if (w(3) >1 && w(3) ~= b(3)) || (w(3) == 1 && length(obj.Bias) ~= 1)
                        error('Inconsistency between filter weights and filter biases');
                    end
                    
                    obj.NumFilters = w(3);
                    obj.NumChannels = w(4);
                    obj.FilterSize = [w(1) w(2)];
                
                case 5
                    
                    layer_name = varargin{1}; 
                    filter_weights = varargin{2};
                    filter_bias = varargin{3};
                    cropping_mat = varargin{4};
                    stride_mat = varargin{5};
                    
                    if ischar(layer_name)
                        obj.Name = layer_name;
                    else
                        error('Layer name should be a char array');
                    end
                    
                    w = size(filter_weights);
                    b = size(filter_bias);
                    
                    if length(w) == 2
                        obj.NumFilters = 1;
                        obj.NumChannels = 1;
                        obj.FilterSize = [w(1) w(2)];
                        obj.Weights = filter_weights;
                        
                     elseif length(w) == 3
                        obj.NumFilters = w(3);
                        obj.NumChannels = 1;
                        obj.FilterSize = [w(1) w(2)];
                        obj.Weights = filter_weights;
                    elseif length(w) == 4
                        obj.NumFilters = w(3);
                        obj.NumChannels = w(4);
                        obj.FilterSize = [w(1) w(2)];
                        obj.Weights = filter_weights;
                    else
                        error('Invalid weight matrix');
                    end
                    
                    if ndims(filter_bias) < 3
                        error('Invalid biases array (must be 3D: [1, 1, NumFilters])');
                    else
                        obj.Bias = filter_bias;
                    end

                    if length(w) == 4 && w(3) ~= b(3)
                        error('Inconsistency between filter weights and filter biases');
                    end

                    if numel(cropping_mat) ~= 4 || ~isrow(cropping_mat)
                        error('Invalid cropping matrix (must be [1x4])');
                    end
                    obj.CroppingSize = cropping_mat;

                    if numel(stride_mat) ~= 2 || ~isrow(stride_mat)
                        error('Invalid stride matrix (must be [1x2])');
                    end
                    obj.Stride = stride_mat;

                case 4

                    filter_weights = varargin{1};
                    filter_bias = varargin{2};
                    cropping_mat = varargin{3};
                    stride_mat = varargin{4};

                    w = size(filter_weights);
                    b = size(filter_bias);

                    if length(w) == 2
                        obj.NumFilters = 1;
                        obj.NumChannels = 1;
                        obj.FilterSize = [w(1) w(2)];
                        obj.Weights = filter_weights;

                    elseif length(w) == 3
                        obj.NumFilters = w(3);
                        obj.NumChannels = 1;
                        obj.FilterSize = [w(1) w(2)];
                        obj.Weights = filter_weights;
                    elseif length(w) == 4
                        obj.NumFilters = w(3);
                        obj.NumChannels = w(4);
                        obj.FilterSize = [w(1) w(2)];
                        obj.Weights = filter_weights;
                    else
                        error('Invalid weight matrix');
                    end

                    % Validate bias: should be [1, 1, NumFilters] but MATLAB drops
                    % trailing singletons, so [1,1,1] becomes [1,1]
                    % Accept 2D [1,1] for single filter case (NumFilters=1)
                    if ndims(filter_bias) == 2 && isequal(size(filter_bias), [1, 1])
                        % Single filter case - 2D [1,1] is acceptable
                        obj.Bias = filter_bias;
                    elseif ndims(filter_bias) >= 3
                        obj.Bias = filter_bias;
                    else
                        error('Invalid biases array (must be 3D: [1, 1, NumFilters])');
                    end

                    if length(w) == 4 && w(3) ~= b(3)
                        error('Inconsistency between filter weights and filter biases');
                    end

                    if numel(cropping_mat) ~= 4 || ~isrow(cropping_mat)
                        error('Invalid cropping matrix (must be [1x4])');
                    end
                    obj.CroppingSize = cropping_mat;

                    if numel(stride_mat) ~= 2 || ~isrow(stride_mat)
                        error('Invalid stride matrix (must be [1x2])');
                    end
                    obj.Stride = stride_mat;

                case 2
                    
                    filter_weights = varargin{1};
                    filter_bias = varargin{2};                
                    
                    w = size(filter_weights);
                    b = size(filter_bias);
                    
                    if length(w) == 2
                        obj.NumFilters = 1;
                        obj.NumChannels = 1;
                        obj.FilterSize = [w(1) w(2)];
                        obj.Weights = filter_weights;
                        
                    elseif length(w) == 3
                        obj.NumFilters = w(3);
                        obj.NumChannels = 1;
                        obj.FilterSize = [w(1) w(2)];
                        obj.Weights = filter_weights;
                    elseif length(w) == 4
                        obj.NumFilters = w(3);
                        obj.NumChannels = w(4);
                        obj.FilterSize = [w(1) w(2)];
                        obj.Weights = filter_weights;
                    else
                        error('Invalid weight matrix');
                    end

                    % Accept 2D [1,1] bias for single filter case
                    % MATLAB drops trailing singletons, so [1,1,1] becomes [1,1]
                    obj.Bias = filter_bias;

                    if length(w) == 4 && w(3) ~= b(3)
                        error('Inconsistency between filter weights and filter biases');
                    end

                otherwise
                    error('Invalid number of inputs (should be 2, 5, or 6)');
                                 
            end

        end
                
    end
    
    
    methods % evaluation methods
                
        % evaluate using matconvnet
        function y = evaluate(obj, input)
            % @input: 3-dimensional array, for example, input(:, :, :)
            % @y: high-dimensional array (output volume), depth of output = number of filters
                        
            %y = vl_nnconvt(input, double(obj.Weights), double(obj.Bias), 'Upsample', obj.Stride, 'Crop', obj.CroppingSize);
            input = dlarray(input);
            y = extractdata(dltranspconv(input,obj.Weights,obj.Bias,"Cropping",obj.CroppingSize,"Stride",obj.Stride,"DataFormat",'SSCU'));

        end
                    
    end
    
       
    methods % reachability analysis 

        % reachability analysis method using ImageStar
        function S = reach_star_single_input(obj, input)
            % @inputs: an ImageStar input set
            % @S: an imagestar with number of channels = obj.NumFilters
            
            if ~isa(input, 'ImageStar')
                error('The input is not an ImageStar object');
            end
            
            if input.numChannel ~= obj.NumChannels
                error("Input set contains %d channels while the convolutional layers has %d channels", input.numChannel, obj.NumChannels);
            end
            
            % compute output sets
            % if isa(input.V, 'gpuArray')
            c = extractdata(dltranspconv(dlarray(input.V(:,:,:,1)), obj.Weights, obj.Bias, "Stride", obj.Stride, "Cropping", obj.CroppingSize,"DataFormat",'SSCU'));
            V = extractdata(dltranspconv(dlarray(input.V(:,:,:,2:input.numPred + 1)), obj.Weights, 0, "Stride", obj.Stride, "Cropping", obj.CroppingSize,"DataFormat",'SSCU'));
            % else
                % c = vl_nnconvt(input.V(:,:,:,1), obj.Weights, obj.Bias, 'Upsample', obj.Stride, 'Crop', obj.CroppingSize);
                % V = vl_nnconvt(input.V(:,:,:,2:input.numPred + 1), obj.Weights, [], 'Upsample', obj.Stride, 'Crop', obj.CroppingSize);
            % end
            Y = cat(4, c, V);
            S = ImageStar(Y, input.C, input.d, input.pred_lb, input.pred_ub);
                  
        end
        
        % reachability analysis using ImageZonos
        function Z = reach_zono(obj, input)
            % @input: an ImageZono input set
            % @Z: an ImageZono with number of channels = obj.NumFilters
            
            if ~isa(input, 'ImageZono')
                error('The input is not an ImageZono object');
            end
            
            if input.numChannels ~= obj.NumChannels
                error("Input set contains %d channels while the convolutional layers has %d channels", input.numChannels, obj.NumChannels);
            end
            
            % compute output sets
            c = vl_nnconvt(input.V(:,:,:,1), obj.Weights, obj.Bias, 'Upsample', obj.Stride, 'Crop', obj.CroppingSize);
            V = vl_nnconvt(input.V(:,:,:,2:input.numPred + 1), obj.Weights, [], 'Upsample', obj.Stride, 'Crop', obj.CroppingSize);         
            Y = cat(4, c, V);
            Z = ImageZono(Y);
            
        end
        
        % multiple inputs (Stars)
        function images = reach_star_multipleInputs(obj, in_images, option)
            % @in_images: an array of ImageStars input set
            % @option: 
            % @images: an array of ImageStars output set
            
            n = length(in_images);
            images(n) = ImageStar; 
            
            if strcmp(option, 'parallel')
                parfor i=1:n
                    images(i) = obj.reach_star_single_input(in_images(i));
                end
            elseif strcmp(option, 'single') || isempty(option)
                for i=1:n
                    images(i) = obj.reach_star_single_input(in_images(i));
                end
            else
                error('Unknown computation option');
            end
            
        end
        
        function images = reach_zono_multipleInputs(obj, in_images, option)
            % @in_images: an array of ImageZonos input set
            % @option: = 'parallel' or 'single' or empty

            n = length(in_images);
            images(n) = ImageZono; 

            if strcmp(option, 'parallel')
                parfor i=1:n
                    images(i) = obj.reach_zono(in_images(i));
                end
            elseif strcmp(option, 'single') || isempty(option)
                for i=1:n
                    images(i) = obj.reach_zono(in_images(i));
                end
            else
                error('Unknown computation option');
            end

        end

        % reach star with multiple inputs
        function images = reach(varargin)
            % @inputs: an array of ImageStar or ImageZono input set
            % @S: an array of ImageStar output set
            
            switch nargin
                case 7
                    obj = varargin{1};
                    in_images = varargin{2};
                    method = varargin{3};
                    option = varargin{4};

                case 6
                    obj = varargin{1};
                    in_images = varargin{2};
                    method = varargin{3};
                    option = varargin{4};

                case 5
                    obj = varargin{1};
                    in_images = varargin{2};
                    method = varargin{3};
                    option = varargin{4};
                case 4
                    obj = varargin{1};
                    in_images = varargin{2};
                    method = varargin{3};
                    option = varargin{4};
                case 3
                    obj = varargin{1};
                    in_images = varargin{2}; 
                    method = varargin{3};
                    option = 'single';
                case 2
                    obj = varargin{1};
                    in_images = varargin{2}; 
                    method = 'approx-star';
                    option = 'single';
                otherwise
                    error('Invalid number of input arguments, should be 1, 2, 3, 4, 5, or 6');
            end
         
            if strcmp(method, 'approx-star') || strcmp(method, 'exact-star') || strcmp(method, 'abs-dom') || contains(method, 'relax-star')
                images = obj.reach_star_multipleInputs(in_images, option);
            elseif strcmp(method, 'approx-zono')
                images = obj.reach_zono_multipleInputs(in_images, option);
            else
                error("Unknown reachability method");
            end
            
        end
        
    end

    methods % helper method

        % change params to gpuArrays
        function obj = toGPU(obj)
            obj.Weights = gpuArray(obj.Weights);
            obj.Bias = gpuArray(obj.Bias);
        end

        % Change params precision
        function obj = changeParamsPrecision(obj, precision)
            if strcmp(precision, "double")
                obj.Weights = double(obj.Weights);
                obj.Bias = double(obj.Bias);
            elseif strcmp(precision, "single")
                obj.Weights = single(obj.Weights);
                obj.Bias = single(obj.Bias);
            else
                error("Parameter numerical precision must be 'single' or 'double'");
            end
        end
        
    end
    
    
    methods(Static) % parsing
        
        % parse a trained transposeConvolutional2dLayer from matlab
        function L = parse(layer)
            % @layer: a transpose convolutional 2d layer from matlab deep
            % neural network tool box
            % @L : a TransposeCov2DLayer for reachability analysis purpose            
            
            if ~isa(layer, 'nnet.cnn.layer.TransposedConvolution2DLayer')
                error('Input is not a Matlab nnet.cnn.layer.TransposedConvolution2DLayer class');
            end
            
            L = TransposedConv2DLayer(layer.Name, layer.Weights, layer.Bias, layer.CroppingSize, layer.Stride, layer.NumInputs, layer.InputNames, layer.NumOutputs, layer.OutputNames);
                                    
        end

    end
    
end

