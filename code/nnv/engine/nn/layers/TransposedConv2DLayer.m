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
        % Cropping = [0 0]; 
        
        
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
            % author: Dung Tran
            % date: 12/5/2018    
            % update: 
            
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
                    if length(b) ~= 3
                        error('Invalid biases array');
                    end

                    if w(3) ~= b(3)
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
                    
                    if length(b) ~= 3
                        error('Invalid biases array');
                    else
                        obj.Bias = filter_bias;
                    end

                    if length(w) == 4 && w(3) ~= b(3)
                        error('Inconsistency between filter weights and filter biases');
                    end
                    
                    p = size(cropping_mat);
                    if length(p) ~= 4 || p(1) ~= 1
                        error('Invalid cropping matrix');
                    end
                    obj.CroppingSize = cropping_mat;
                                               
                    s = size(stride_mat);
                    if length(s) ~= 2 || s(1) ~= 1
                        error('Invalid stride matrix');
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
                    
                    if length(b) ~= 3
                        error('Invalid biases array');
                    else
                        obj.Bias = filter_bias;
                    end

                    if length(w) == 4 && w(3) ~= b(3)
                        error('Inconsistency between filter weights and filter biases');
                    end
                    
                    p = size(cropping_mat);
                    if length(p) ~= 4 || p(1) ~= 1
                        error('Invalid cropping matrix');
                    end
                    obj.CroppingSize = cropping_mat;
                                               
                    s = size(stride_mat);
                    if length(s) ~= 2 || s(1) ~= 1
                        error('Invalid stride matrix');
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

                    obj.Bias = filter_bias;

                    if length(w) == 4 && w(3) ~= b(3)
                        error('Inconsistency between filter weights and filter biases');
                    end
                                        
                           
                                    
                otherwise
                    error('Invalid number of inputs (should be 2, 5, or 6)');
                                 
            end
                    

        end
        
                
    end
    
    
    % evaluation method
    methods
                
        % evaluate using matconvnet
        function y = evaluate(obj, input)
            % @input: 3-dimensional array, for example, input(:, :, :)
            % @y: high-dimensional array (output volume), depth of output = number of filters
            % @option: 'single' or 'double' or empty precision of
            % computation
            
            
            % author: Dung Tran
            % date: 4/22/2020
                        
            y = vl_nnconvt(input, double(obj.Weights), double(obj.Bias), 'Upsample', obj.Stride, 'Crop', obj.CroppingSize);

        end
                    

        
    end
    
       
    % reachability analysis using star set
    
    methods
        % reachability analysis method using ImageStar
        function S = reach_star_single_input(obj, input)
            % @inputs: an ImageStar input set
            % @S: an imagestar with number of channels = obj.NumFilters
            
            % author: Dung Tran
            % date: 6/11/2019
            
            if ~isa(input, 'ImageStar')
                error('The input is not an ImageStar object');
            end
            
            if input.numChannel ~= obj.NumChannels
                error("Input set contains %d channels while the convolutional layers has %d channels", input.numChannel, obj.NumChannels);
            end
            
            % compute output sets
            c = vl_nnconvt(double(input.V(:,:,:,1)), double(obj.Weights), double(obj.Bias), 'Upsample', obj.Stride, 'Crop', obj.CroppingSize);
            V = vl_nnconvt(double(input.V(:,:,:,2:input.numPred + 1)), double(obj.Weights), [], 'Upsample', obj.Stride, 'Crop', obj.CroppingSize);         
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
            c = vl_nnconvt(double(input.V(:,:,:,1)), double(obj.Weights), double(obj.Bias), 'Upsample', obj.Stride, 'Crop', obj.CroppingSize);
            V = vl_nnconvt(double(input.V(:,:,:,2:input.numPred + 1)), double(obj.Weights), [], 'Upsample', obj.Stride, 'Crop', obj.CroppingSize);         
            Y = cat(4, c, V);
            Z = ImageZono(Y);
            
        end
        
        
        % hangle multiple inputs
        function images = reach_star_multipleInputs(obj, in_images, option)
            % @in_images: an array of ImageStars input set
            % @option: 
            % @images: an array of ImageStars output set
            
            % author: Dung Tran
            % date: 1/7/2020
            
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

            % author: Dung Tran
            % date: 1/7/2020

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
                        
            % author: Dung Tran
            % date: 7/16/2019
            
            switch nargin
                case 7
                    obj = varargin{1};
                    in_images = varargin{2};
                    method = varargin{3};
                    option = varargin{4};
                    % relaxFactor = varargin{5}; do not use
                    % dis_opt = varargin{6}; do not use
                    % lp_solver = varargin{7}; do not use
                case 6
                    obj = varargin{1};
                    in_images = varargin{2};
                    method = varargin{3};
                    option = varargin{4};
                    %relaxFactor = varargin{5}; do not use
                    % dis_opt = varargin{6}; do not use
                case 5
                    obj = varargin{1};
                    in_images = varargin{2};
                    method = varargin{3};
                    option = varargin{4};
                    %relaxFactor = varargin{5}; do not use
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
    
    
    methods(Static) % parsing, get zero input pading, compute feature maps
        
        % parse a trained transposeConvolutional2dLayer from matlab
        function L = parse(layer)
            % @layer: a transpose convolutional 2d layer from matlab deep
            % neural network tool box
            % @L : a TransposeCov2DLayer for reachability analysis purpose
            
            % author: Dung Tran
            % date: 4/22/2020
            
            
            if ~isa(layer, 'nnet.cnn.layer.TransposedConvolution2DLayer')
                error('Input is not a Matlab nnet.cnn.layer.TransposedConvolution2DLayer class');
            end
            
            L = TransposedConv2DLayer(layer.Name, layer.Weights, layer.Bias, layer.CroppingSize, layer.Stride, layer.NumInputs, layer.InputNames, layer.NumOutputs, layer.OutputNames);
                        
            fprintf('\nParsing a Matlab transposed convolutional 2d layer is done successfully');
            
        end

    end
    
end

