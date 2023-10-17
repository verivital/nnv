classdef Conv3DLayer < handle
    % The convolutional 3D layer class
    %   Contain constructor and reachability analysis methods
    % Main references:
    %    https://www.mathworks.com/help/deeplearning/ref/nnet.cnn.layer.convolution3dlayer.html
    %
    %   Diego Manzanas: 10/11/2023
    
    properties
        Name = 'convolutional_3d_layer';
        % Hyperparameters
        FilterSize = []; % height and width of filters
        NumChannels = 'auto';
        NumFilters = 0; % number of filters
        Stride = [1 1 1]; % step size for traversing input
        DilationFactor = [1 1 1]; % factor for dilated convolution
        PaddingMode = 'manual';
        PaddingSize = [0 0 0 0]; % size of padding [t b l r] for nonnegative integers
        
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
        function obj = Conv3DLayer(varargin)
            
            switch nargin
                
                case 10 % used for parsing a Matlab conv3dlayer 
                    
                    layer_name = varargin{1}; 
                    filter_weights = varargin{2};
                    filter_bias = varargin{3};
                    padding_mat = varargin{4};
                    stride_mat = varargin{5};
                    dilation_mat = varargin{6};
                    obj.NumInputs = varargin{7};
                    obj.InputNames = varargin{8};
                    obj.NumOutputs = varargin{9};
                    obj.OutputNames = varargin{10};
                    
                    if ischar(layer_name)
                        obj.Name = layer_name;
                    else
                        error('Layer name should be a char array');
                    end
                    
                    w = size(filter_weights);
                    b = size(filter_bias);

                    if length(w) == 3
                        obj.NumFilters = 1;
                        obj.NumChannels = 1;
                        obj.FilterSize = [w(1) w(2) w(3)];
                        obj.Weights = filter_weights;
                        
                    elseif length(w) == 4
                        obj.NumFilters = 1;
                        obj.NumChannels = w(4);
                        obj.FilterSize = [w(1) w(2) w(3)];
                        obj.Weights = filter_weights;

                    elseif length(w) == 5
                        obj.NumFilters = w(5);
                        obj.NumChannels = w(4);
                        obj.FilterSize = [w(1) w(2) w(3)];
                        obj.Weights = filter_weights;

                    else
                        error('Invalid weight matrix');
                    end

                    if length(w) == 5 && w(5) ~= b(4)
                        error('Inconsistency between filter weights and filter biases');
                    end

                    obj.Bias = filter_bias;
                                        
                    p = size(padding_mat);
                   
                    if p(2) ~= 3|| p(1) ~= 2
                        error('Invalid padding matrix');
                    end
                    obj.PaddingSize = padding_mat;
                    
                    s = size(stride_mat);
                    if s(2) ~= 3 || s(1) ~= 1
                        error('Invalid stride matrix');
                    end
                    obj.Stride = stride_mat;
                    
                    d = size(dilation_mat);
                    if d(2) ~= 3 || d(1) ~= 1
                        error('Invalid dilation matrix');
                    end                
                    obj.DilationFactor = dilation_mat;
                    
                case 6
                    
                    layer_name = varargin{1}; 
                    filter_weights = varargin{2};
                    filter_bias = varargin{3};
                    padding_mat = varargin{4};
                    stride_mat = varargin{5};
                    dilation_mat = varargin{6};
                    
                    if ischar(layer_name)
                        obj.Name = layer_name;
                    else
                        error('Layer name should be a char array');
                    end
                    
                    w = size(filter_weights);
                    b = size(filter_bias);

                    if length(w) == 3
                        obj.NumFilters = 1;
                        obj.NumChannels = 1;
                        obj.FilterSize = [w(1) w(2) w(3)];
                        obj.Weights = filter_weights;
                        
                    elseif length(w) == 4
                        obj.NumFilters = 1;
                        obj.NumChannels = w(4);
                        obj.FilterSize = [w(1) w(2) w(3)];
                        obj.Weights = filter_weights;

                    elseif length(w) == 5
                        obj.NumFilters = w(5);
                        obj.NumChannels = w(4);
                        obj.FilterSize = [w(1) w(2) w(3)];
                        obj.Weights = filter_weights;

                    else
                        error('Invalid weight matrix');
                    end
                    
                    if length(w) == 5 && w(5) ~= b(4)
                        error('Inconsistency between filter weights and filter biases');
                    end

                    obj.Bias = filter_bias;
                                        
                    p = size(padding_mat);
                   
                    if p(2) ~= 3|| p(1) ~= 2
                        error('Invalid padding matrix');
                    end
                    obj.PaddingSize = padding_mat;
                    
                    s = size(stride_mat);
                    if s(2) ~= 3 || s(1) ~= 1
                        error('Invalid stride matrix');
                    end
                    obj.Stride = stride_mat;
                    
                    d = size(dilation_mat);
                    if d(2) ~= 2 || d(1) ~= 1
                        error('Invalid dilation matrix');
                    end                
                    obj.DilationFactor = dilation_mat;
                    
                case 5
                    
                    filter_weights = varargin{1};
                    filter_bias = varargin{2};
                    padding_mat = varargin{3};
                    stride_mat = varargin{4};
                    dilation_mat = varargin{5};
                    
                    w = size(filter_weights);
                    b = size(filter_bias);

                    if length(w) == 3
                        obj.NumFilters = 1;
                        obj.NumChannels = 1;
                        obj.FilterSize = [w(1) w(2) w(3)];
                        obj.Weights = filter_weights;
                        
                    elseif length(w) == 4
                        obj.NumFilters = 1;
                        obj.NumChannels = w(4);
                        obj.FilterSize = [w(1) w(2) w(3)];
                        obj.Weights = filter_weights;

                    elseif length(w) == 5
                        obj.NumFilters = w(5);
                        obj.NumChannels = w(4);
                        obj.FilterSize = [w(1) w(2) w(3)];
                        obj.Weights = filter_weights;

                    else
                        error('Invalid weight matrix');
                    end
                    
                    if length(w) == 5 && w(5) ~= b(4)
                        error('Inconsistency between filter weights and filter biases');
                    end

                    obj.Bias = filter_bias;
                                        
                    p = size(padding_mat);
                   
                    if p(2) ~= 3|| p(1) ~= 2
                        error('Invalid padding matrix');
                    end
                    obj.PaddingSize = padding_mat;
                    
                    s = size(stride_mat);
                    if s(2) ~= 3 || s(1) ~= 1
                        error('Invalid stride matrix');
                    end
                    obj.Stride = stride_mat;
                    
                    d = size(dilation_mat);
                    if d(2) ~= 2 || d(1) ~= 1
                        error('Invalid dilation matrix');
                    end                
                    obj.DilationFactor = dilation_mat;                 
                    
                case 2
                    
                    filter_weights = varargin{1};
                    filter_bias = varargin{2};
                    
                    w = size(filter_weights);
                    b = size(filter_bias);
                    
                    if length(w) == 3
                        obj.NumFilters = 1;
                        obj.NumChannels = 1;
                        obj.FilterSize = [w(1) w(2) w(3)];
                        obj.Weights = filter_weights;
                        
                    elseif length(w) == 4
                        obj.NumFilters = 1;
                        obj.NumChannels = w(4);
                        obj.FilterSize = [w(1) w(2) w(3)];
                        obj.Weights = filter_weights;

                    elseif length(w) == 5
                        obj.NumFilters = w(5);
                        obj.NumChannels = w(4);
                        obj.FilterSize = [w(1) w(2) w(3)];
                        obj.Weights = filter_weights;

                    else
                        error('Invalid weight matrix');
                    end           

                    obj.Bias = filter_bias;

                    if length(w) == 5 && w(5) ~= b(4)
                        error('Inconsistency between filter weights and filter biases');
                    end
                                        
                    % use default setting for Stride, Padding and Dilation
                                    
                otherwise
                    error('Invalid number of inputs (should be 2, 5, or 6)');
                                 
            end
             
        end
        
        % set stride method 
        function set_stride(obj, stride)
            % @stride: stride matrix
            
            [n, m] = size(stride);
            
            if n ~= 1
                error('Stride matrix shoule have one row');
            end
            
            if m == 1
                obj.Stride = [stride stride stride];
            elseif m == 3
                obj.Stride = stride;
            elseif m ~= 1 && m ~= 3
                error('Invalide stride matrix');
            end
            
        end
        
        % set dilation factor
        function set_dilation(obj, dilation)
            % @dilation: dilation matrix
            
            [n, m] = size(dilation);
            
            if n ~= 1
                error('Dilation matrix shoule have one row');
            end
            
            if m == 1
                obj.DilationFactor = [dilation dilation dilation];
            elseif m == 3
                obj.DilationFactor = dilation;
            elseif m ~= 1 && m ~= 3
                error('Invalide dilation matrix');
            end
            
        end
        
        % set padding 
        function set_padding(obj, padding)
            % @padding: padding matrix
            
            [n, m] = size(padding);
            
            if n ~= 2
                error('Padding matrix shoule have two rows.');
            end
            
            if m == 2
                obj.PaddingSize = padding * ones(2,3);
            elseif m == 3
                obj.PaddingSize = padding;
            elseif m ~= 2 && m ~= 3
                error('Invalide padding matrix');
            end
        end
        
        % set weights and biases
        function set_weights_biases(obj, weights, biases)
            % @weights: a 5-dimensional array
            %    [FilterSize(1)-by-FilterSize(2)-by-FilterSize(3)-by-NumChannels-by-NumFilter]
            % @biases: a 4-dimensional array: [1-by-1-by-1-by-NumFilters]
            
            w = size(weights);
            
            if length(w) == 3
                obj.NumFilters = 1;
                obj.NumChannels = 1;
            elseif length(w) == 4
                obj.NumFilters = 1;
                obj.NumChannels = w(4);
            elseif length(w) == 5
                obj.NumFilters = w(5);
                obj.NumChannels = w(4);
            else
                error('Invalid weights array');
            end
            
            obj.Weights = weights;
            obj.Bias = biases;
            
        end
        
        % set name of this layer
        function set_name(obj, name)
            % @name: name of this layer
            
            if ischar(name)
                obj.Name = name;
            else
                error('Input is not a character array');
            end
            
        end
        
    end
    
    
    methods % evaluation methods
        
        % evaluate (using dlconv)
        function y = evaluate(obj, input) 
            % @input: 4-dimensional array, for example, input(:, :, :, :)
            % @y: high-dimensional array (output volume), depth of output = number of filters
            
            x = dlarray(input, "SSSC");
            y = dlconv(x, obj.Weights, obj.Bias, 'Stride', obj.Stride, 'Padding', obj.PaddingSize, 'DilationFactor', obj.DilationFactor);
            y = extractdata(y);

        end

    end
    
    
    methods % reachability analysis functions

        % an VolumeStar represent a set of 3d images (volumes) (3D matrix of h x w x d)
        function S = reach_star_single_input(obj, input) 
            % @inputs: a VolumeStar input set
            % @S: a VolumeStar with number of channels = obj.NumFilters
            
            if ~isa(input, 'VolumeStar')
                error('The input is not an VolumeStar object');
            end
            
            if input.numChannel ~= obj.NumChannels
                error("Input set contains %d channels while the convolutional layers has %d channels", input.numChannel, obj.NumChannels);
            end
            
            % compute output sets
            c = dlconv(dlarray(input.V(:,:,:,:,1), "SSSC"), obj.Weights, obj.Bias, 'Stride', obj.Stride, 'Padding', obj.PaddingSize, 'DilationFactor', obj.DilationFactor);
            V = dlconv(dlarray(input.V(:,:,:,:,2:input.numPred + 1), "SSSCB"), obj.Weights, 0, 'Stride', obj.Stride, 'Padding', obj.PaddingSize, 'DilationFactor', obj.DilationFactor);        
            Y = cat(5, extractdata(c), extractdata(V));
            S = VolumeStar(Y, input.C, input.d, input.pred_lb, input.pred_ub);
                  
        end
        
        % handle multiple inputs (VolumeStars)
        function volumes = reach_star_multipleInputs(obj, in_volumes, option)
            % @in_ivolumes: an array of VolumeStars input set
            % @option: {'single', 'parallel'} 
            % @volumes: an array of VolumeStars output set
            
            n = length(in_volumes);
            volumes(n) = VolumeStar; 
            
            if strcmp(option, 'parallel')
                parfor i=1:n
                    volumes(i) = obj.reach_star_single_input(in_volumes(i));
                end
            elseif strcmp(option, 'single') || isempty(option)
                for i=1:n
                    volumes(i) = obj.reach_star_single_input(in_volumes(i));
                end
            else
                error('Unknown computation option');

            end
            
        end
        
        % Main reach function
        function volumes = reach(varargin)
            % @inputs: an array of VolumeStar
            % @volumes: an array of VolumeStar output set
            
            switch nargin
                
                 case 7
                    obj = varargin{1};
                    in_volumes = varargin{2};
                    method = varargin{3};
                    option = varargin{4};
                    % relaxFactor = varargin{5}; do not use
                    % dis_opt = varargin{6}; do not use
                    % lp_solver = varargin{7}; do not use
                case 6
                    obj = varargin{1};
                    in_volumes = varargin{2};
                    method = varargin{3};
                    option = varargin{4};
                    %relaxFactor = varargin{5}; do not use
                    % dis_opt = varargin{6}; do not use
                
                case 5
                    obj = varargin{1};
                    in_volumes = varargin{2};
                    method = varargin{3};
                    option = varargin{4};
                    %relaxFactor = varargin{5}; do not use
                
                case 4
                    obj = varargin{1};
                    in_volumes = varargin{2};
                    method = varargin{3};
                    option = varargin{4};
                case 3
                    obj = varargin{1};
                    in_volumes = varargin{2}; 
                    method = varargin{3};
                    option = 'single';
                case 2
                    obj = varargin{1};
                    in_volumes = varargin{2}; 
                    method = 'approx-star';
                    option = 'single';
                    
                otherwise
                    error('Invalid number of input arguments, should be 1, 2, 3, 4, 5 or 6');
            end
         
            if strcmp(method, 'approx-star') || strcmp(method, 'exact-star') || strcmp(method, 'abs-dom')|| contains(method, "relax-star")
                volumes = obj.reach_star_multipleInputs(in_volumes, option);              
            else
                error("Unknown reachability method");
            end
            
        end

    end

    methods(Static) % parsing
        
        % parse a trained convolutional3dLayer from matlab
        function L = parse(conv3dLayer)
            % @conv3dLayer: a convolutional 3d layer from matlab deep
            % neural network tool box
            % @L : a Conv3DLayer for reachability analysis purpose

            if ~isa(conv3dLayer, 'nnet.cnn.layer.Convolution3DLayer')
                error('Input is not a Matlab nnet.cnn.layer.Convolution3DLayer class');
            end
            
            layer_name = conv3dLayer.Name; 
            filter_weights = conv3dLayer.Weights;
            filter_bias = conv3dLayer.Bias;
            padding_mat = conv3dLayer.PaddingSize;
            stride_mat = conv3dLayer.Stride;
            dilation_mat = conv3dLayer.DilationFactor;
            
            L = Conv3DLayer(layer_name, filter_weights, filter_bias, padding_mat, stride_mat, dilation_mat, conv3dLayer.NumInputs, conv3dLayer.InputNames, conv3dLayer.NumOutputs, conv3dLayer.OutputNames);
                                    
        end
        
    end
    
end

