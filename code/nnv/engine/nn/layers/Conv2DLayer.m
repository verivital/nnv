classdef Conv2DLayer < handle
    % The convolutional 2D layer class
    %   Contain constructor and reachability analysis methods
    % Main references:
    % 1) An intuitive explanation of convolutional neural networks: 
    %    https://ujjwalkarn.me/2016/08/11/intuitive-explanation-convnets/
    % 2) More detail about mathematical background of CNN
    %    http://cs231n.github.io/convolutional-networks/
    % 3) Matlab implementation of Convolution2DLayer (for training and evaluating purpose)
    %    https://www.mathworks.com/help/deeplearning/ug/layers-of-a-convolutional-neural-network.html
    %    https://www.mathworks.com/help/deeplearning/ref/nnet.cnn.layer.convolution2dlayer.html
    
    %   Dung Tran: 12/5/2018
    
    properties
        Name = 'convolutional_layer';
        % Hyperparameters
        FilterSize = []; % height and width of filters
        NumChannels = 'auto';
        NumFilters = 0; % number of filters
        Stride = [1 1]; % step size for traversing input
        DilationFactor = [1 1]; % factor for dilated convolution
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
        function obj = Conv2DLayer(varargin)           
            % author: Dung Tran
            % date: 12/5/2018    
            % update: 
            
            switch nargin
                
                case 10 % used for parsing a Matlab conv2dlayer 
                    
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
                    
                    if length(w) ~= 4
                        error('Invalid weights array');
                    end
                    if length(b) ~= 3
                        error('Invalid biases array');
                    end

                    if w(4) ~= b(3)
                        error('Inconsistency between filter weights and filter biases');
                    end
                    
                    obj.NumFilters = w(4);
                    obj.NumChannels = w(3);
                    obj.FilterSize = [w(1) w(2)];
                    obj.Weights = filter_weights;
                    obj.Bias = filter_bias;
                                        
                    p = size(padding_mat);
                   
                    if  length(p) ~= 2 || p(2) ~= 4|| p(1) ~= 1
                        error('Invalid padding matrix');
                    end
                    obj.PaddingSize = padding_mat;
                    
                    s = size(stride_mat);
                    if length(s) ~= 2 || s(1) ~= 1
                        error('Invalid stride matrix');
                    end
                    obj.Stride = stride_mat;
                    
                    d = size(dilation_mat);
                    if length(d) ~= 2 || d(1) ~= 1
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
                    
                    if length(w) ~= 4
                        error('Invalid weights array');
                    end
                    if length(b) ~= 3
                        error('Invalid biases array');
                    end

                    if w(4) ~= b(3)
                        error('Inconsistency between filter weights and filter biases');
                    end
                    
                    obj.NumFilters = w(4);
                    obj.NumChannels = w(3);
                    obj.FilterSize = [w(1) w(2)];
                    obj.Weights = filter_weights;
                    obj.Bias = filter_bias;
                                        
                    p = size(padding_mat);
                   
                    if  length(p) ~= 2 || p(2) ~= 4|| p(1) ~= 1
                        error('Invalid padding matrix');
                    end
                    obj.PaddingSize = padding_mat;
                    
                    s = size(stride_mat);
                    if length(s) ~= 2 || s(1) ~= 1
                        error('Invalid stride matrix');
                    end
                    obj.Stride = stride_mat;
                    
                    d = size(dilation_mat);
                    if length(d) ~= 2 || d(1) ~= 1
                        error('Invalid dilation matrix');
                    end                
                    obj.DilationFactor = dilation_mat;
                    
                    
                
                case 5
                    
                    filter_weights = varargin{1};
                    filter_bias = varargin{2};
                    padding_mat = varargin{3};
                    stride_mat = varargin{4};
                    dilation_mat = varargin{5};
                    
                    obj.Name = 'convolutional_layer';                   
                    
                    w = size(filter_weights);
                    b = size(filter_bias);
                    
                    if length(w) == 2
                        obj.NumFilters = 1;
                        obj.NumChannels = 1;
                        obj.FilterSize = [w(1) w(2)];
                        obj.Weights = filter_weights;
                        
                    elseif length(w) == 3
                        obj.NumFilters = 1;
                        obj.NumChannels = w(3);
                        obj.FilterSize = [w(1) w(2)];
                        obj.Weights = filter_weights;
                    elseif length(w) == 4
                        obj.NumFilters = w(4);
                        obj.NumChannels = w(3);
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

                    if length(w) == 4 && w(4) ~= b(3)
                        error('Inconsistency between filter weights and filter biases');
                    end
                    
                    p = size(padding_mat);
                    if length(p) ~= 4 || p(1) ~= 1
                        error('Invalid padding matrix');
                    end
                    obj.PaddingSize = padding_mat;
                                               
                    s = size(stride_mat);
                    if length(s) ~= 2 || s(1) ~= 1
                        error('Invalid stride matrix');
                    end
                    obj.Stride = stride_mat;
                    
                    d = size(dilation_mat);
                    if length(d) ~= 2 || d(1) ~= 1
                        error('Invalid dilation matrix');
                    end                
                    obj.DilationFactor = dilation_mat;                   
                    
                    
                case 2
                    
                    filter_weights = varargin{1};
                    filter_bias = varargin{2};
                                        
                    obj.Name = 'convolutional_layer';                   
                    
                    w = size(filter_weights);
                    b = size(filter_bias);
                    
                    if length(w) == 2
                        obj.NumFilters = 1;
                        obj.NumChannels = 1;
                        obj.FilterSize = [w(1) w(2)];
                        obj.Weights = filter_weights;
                        
                    elseif length(w) == 3
                        obj.NumFilters = 1;
                        obj.NumChannels = w(3);
                        obj.FilterSize = [w(1) w(2)];
                        obj.Weights = filter_weights;
                    elseif length(w) == 4
                        obj.NumFilters = w(4);
                        obj.NumChannels = w(3);
                        obj.FilterSize = [w(1) w(2)];
                        obj.Weights = filter_weights;
                    else
                        error('Invalid weight matrix');
                    end                  

                    obj.Bias = filter_bias;

                    if length(w) == 4 && w(4) ~= b(3)
                        error('Inconsistency between filter weights and filter biases');
                    end
                                        
                    % use default setting for Stride, Padding and Dilation
                    obj.Stride = [1 1]; % step size for traversing input
                    obj.DilationFactor = [1 1]; % factor for dilated convolution
                    obj.PaddingMode = 'manual';
                    obj.PaddingSize = [0 0 0 0]; % size of padding [t b l r] for nonnegative integers
        
                                    
                otherwise
                    error('Invalid number of inputs (should be 2, 5, or 6)');
                                 
            end
                    
                
            
            
            
            
             
        end
        
        % set stride method 
        function set_stride(obj, stride)
            % @stride: stride matrix
            % author: Dung Tran
            % date: 12/5/2018
            
            [n, m] = size(stride);
            
            if n ~= 1
                error('Stride matrix shoule have one row');
            end
            
            if m == 1
                obj.Stride = [stride stride];
            elseif m == 2
                obj.Stride = [stride(1) stride(2)];
            elseif m ~=1 && m ~=2 
                error('Invalide stride matrix');
            end
            
        end
        
        
        % set dilation factor
        function set_dilation(obj, dilation)
            % @dilation: dilation matrix
            % author: Dung Tran
            % date: 12/5/2018
            
            [n, m] = size(dilation);
            
            if n ~= 1
                error('Dilation matrix shoule have one row');
            end
            
            if m == 1
                obj.DilationFactor = [dilation dilation];
            elseif m == 2
                obj.DilationFactor = [dilation(1) dilation(2)];
            elseif m ~=1 && m ~=2 
                error('Invalide dilation matrix');
            end
            
        end
        
        % set padding 
        function set_padding(obj, padding)
            % @padding: padding matrix
            % author: Dung Tran
            % date: 12/5/2018
            
            [n, m] = size(padding);
            
            if n ~= 1
                error('Padding matrix shoule have one row');
            end
            
            if m == 1
                obj.PaddingSize = [padding padding padding padding];
            elseif m == 4
                obj.PaddingSize = [padding(1) padding(2) padding(3) padding(4)];
            elseif m ~=1 && m ~=4 
                error('Invalide padding matrix');
            end
        end
        
        % set weights and biases
        function set_weights_biases(obj, weights, biases)
            % @weights: a 4-dimensional array
            %  = FilterSize(1)-by-FilterSize(2)-by-NumChannels-by-NumFilter
            % @biases: a 3-dimensional array
            %  = 1-by-1-by-NumFilters 
            
            w = size(weights);
            b = size(biases);
            
            if length(w) == 2
                obj.NumFilters = 1;
                obj.NumChannels = 1;
            elseif length(w) == 3
                obj.NumFilters = 1;
                obj.NumChannels = w(3);
            elseif length(w) == 4
                obj.NumFilters = w(4);
                obj.NumChannels = w(3);
            else
                error('Invalid weights array');
            end
            
            obj.Weights = weights;
            obj.Bias = biases;
            
        end
        
        % set name of this layer
        function set_name(obj, name)
            % @name: name of this layer
            
            % author: Dung Tran
            % date: 12/5/2018
            
            if ischar(name)
                obj.Name = name;
            else
                error('Input is not a character array');
            end
            
        end
        
    end
    
    
    % evaluation method
    methods
        
        function y = evaluate2(obj, input, option)
            % @input: 3-dimensional array, for example, input(:, :, :)
            % @y: high-dimensional array (output volume), depth of output = number of filters
            % @option = 'single'
            
            % author: Dung Tran
            % date: 12/10/2018
            
            I1 = input(:,:,1);
            W1 = obj.Weights(:,:,1,1);
            if strcmp(option, 'single') && ~strcmp(option, 'double')
                obj.Weights = single(obj.Weights);
                obj.Bias = single(obj.Bias);
            elseif strcmp(option, 'double')
                obj.Weights = double(obj.Weights);
                obj.Bias = double(obj.Bias);
            else
                error('Unknown precison option');               
            end
            [h, w] = Conv2DLayer.get_size_featureMap(I1, W1, obj.PaddingSize, obj.Stride, obj.DilationFactor);   
            y(:,:, obj.NumFilters) = zeros(h, w, option); % preallocate 3D matrix
            for i=1:obj.NumFilters % number of filters                   
                y(:, :, i) = obj.Bias(:, :, i)* ones(h, w, option); % initialize output by bias matrix    
                % compute feature map with i^th filter 
                for j=1:obj.NumChannels % filter for j^th channel of input 
                    W1 = obj.Weights(:,:,j, i);   
                    I1 = input(:, :, j); % the j^th input channel
                    if j==1
                        y(:, :, i) = Conv2DLayer.compute_featureMap(I1, W1, obj.PaddingSize, obj.Stride, obj.DilationFactor);
                    else
                        y(:, :, i) = y(:, :, i) + Conv2DLayer.compute_featureMap(I1, W1, obj.PaddingSize, obj.Stride, obj.DilationFactor);
                    end                    
                end
                [ny, my] = size(y(:, :, i));
                y(:, :, i) = y(:, :, i) + double(obj.Bias(:, :, i)) * ones(ny, my, option);                               
            end
                   
        end
        
        % evaluate using matconvnet
        function y = evaluate(varargin)
            % @input: 3-dimensional array, for example, input(:, :, :)
            % @y: high-dimensional array (output volume), depth of output = number of filters
            % @option: 'single' or 'double' or empty precision of
            % computation
            
            
            % author: Dung Tran
            % date: 7/18/2019
            
            switch nargin
                
                case 2
                    obj = varargin{1};
                    input = varargin{2};
                    option = 'double';
                    input = double(input);

                case 3
                    obj = varargin{1};
                    input = varargin{2};
                    option = varargin{3};
                                      
                otherwise 
                    error('Invalid number of inputs, should be 1 or 2');
                
            end
            
            if strcmp(option, 'single')
                obj.Weights = single(obj.Weights);
                obj.Bias = single(obj.Bias);
            elseif strcmp(option, 'double')
                obj.Weights = double(obj.Weights);
                obj.Bias = double(obj.Bias);
            else
                error('Unknown precision option');
            end
            
            y = vl_nnconv(input, obj.Weights, obj.Bias, 'Stride', obj.Stride, 'Pad', obj.PaddingSize, 'Dilate', obj.DilationFactor);

        end
                    

        
    end
    
    
    
    methods(Static) % parsing, get zero input pading, compute feature maps
        
        % parse a trained convolutional2dLayer from matlab
        function L = parse(conv2dLayer)
            % @conv2dLayer: a convolutional 2d layer from matlab deep
            % neural network tool box
            % @L : a Cov2DLayer for reachability analysis purpose
            
            % author: Dung Tran
            % date: 12/5/2018
            
            
            if ~isa(conv2dLayer, 'nnet.cnn.layer.Convolution2DLayer')
                error('Input is not a Matlab nnet.cnn.layer.Convolution2DLayer class');
            end
            
            
            layer_name = conv2dLayer.Name; 
            filter_weights = conv2dLayer.Weights;
            filter_bias = conv2dLayer.Bias;
            padding_mat = conv2dLayer.PaddingSize;
            stride_mat = conv2dLayer.Stride;
            dilation_mat = conv2dLayer.DilationFactor;
            
            L = Conv2DLayer(layer_name, filter_weights, filter_bias, padding_mat, stride_mat, dilation_mat, conv2dLayer.NumInputs, conv2dLayer.InputNames, conv2dLayer.NumOutputs, conv2dLayer.OutputNames);
                        
            fprintf('\nParsing a Matlab convolutional 2d layer is done successfully');
            
        end
        
        % parse input image with padding
        function I = get_zero_padding_input(input, paddingSize)
            % @input: an array of input image, 1 or high-dimensional array
            % @paddingSize: paddingSize to construct the new input I
            % @I: the new input array affer applying padding
            
            % author: Dung Tran
            % date: 12/10/2018
            
            n = size(input);        
            if length(paddingSize) ~= 4
                error("\Invalid padding size");
            end
            
            t = paddingSize(1);
            b = paddingSize(2);
            l = paddingSize(3);
            r = paddingSize(4);
 
            
            if length(n) == 2 
                % input volume has only one channel                
                h = n(1); % height of input
                w = n(2); % width of input 
                      
                I = zeros(t + h + b, l + w + r);
                I(t+1:t+h,l+1:l+w) = input; % new constructed input volume
                
            elseif length(n) > 2
                % input volume may has more than one channel
                            
                h = n(1); % height of input
                w = n(2); % width of input 
                d = n(3); % depth of input (number of channel)
               
                I = zeros(t + h + b, l + w + r, d); % preallocate 3D matrix
                for i=1:d
                    I(t+1:t+h, l+1:l+w, i) = input(:, :, i); % new constructed input volume
                end              
                
                
            else
                error('Invalid input');
            end
                         
        end
        
        
        % compute feature map for specific input and weight
        function featureMap = compute_featureMap(input, W, padding, stride, dilation)
            % @input: is input 
            % @W: is a weight matrix (filter)
            % @padding: zero-padding size
            % @stride: step size for traversing input
            % @dilation: factor for dilated convolution     
            % @featureMap: convolved feature (also called feature map)
            
            
            % author: Dung Tran
            % date: 12/10/2018
            
            % referece: 1) https://ujjwalkarn.me/2016/08/11/intuitive-explanation-convnets/
            %           2) https://www.mathworks.com/help/deeplearning/ug/layers-of-a-convolutional-neural-network.html
            
            
            I = Conv2DLayer.get_zero_padding_input(input, padding);
            n = size(I); % n(1) is height and n(2) is width of input
            m = size(W); % m(1) is height and m(2) is width of the filter
            
            % I, W is 2D matrix
            % I is assumed to be the input after zero padding
            % output size: 
            % (InputSize - (FilterSize - 1)*Dilation + 1)/Stride
            
            h = floor((n(1) - (m(1) - 1) * dilation(1) - 1) / stride(1) + 1);  % height of feature map
            w = floor((n(2) - (m(2) - 1) * dilation(2) - 1) / stride(2) + 1);  % width of feature map

            % a collection of start points (the top-left corner of a square) of the region of input that is filtered
            map = cell(h, w); 
            
            for i=1:h
                for j=1:w
                    map{i, j} = zeros(1, 2);

                    if i==1
                        map{i, j}(1) = 1;
                    end
                    if j==1
                        map{i, j}(2) = 1;
                    end

                    if i > 1
                        map{i, j}(1) = map{i - 1, j}(1) + stride(1);
                    end

                    if j > 1
                        map{i, j}(2) = map{i, j - 1}(2) + stride(2);
                    end

                end
            end
            
            % compute feature map for each cell of map
            % do it in parallel using cpu or gpu
            % TODO: explore power of using GPU for this problem
            
            featureMap = zeros(1, h*w); % using single vector for parallel computation
           

            for l=1:h*w
                a = mod(l, w);
                if a == 0
                    i = floor(l/w);
                    j = w;
                else
                    j = a;
                    i = floor(l/w) + 1;
                end

                % get a filtered region
                val = 0;
                i0 = map{i, j}(1);
                j0 = map{i, j}(2);

                ie = i0;
                je = j0;
                for i=1:m(1)
                    for j=1:m(2)
                        if ie <= n(1) && je <= n(2)
                            val = val + I(ie, je) * W(i, j);
                        end
                        je = je + dilation(2);
                    end
                    je = j0;
                    ie = ie + dilation(1);
                end

                featureMap(1,l) = val;

            end

            featureMap = reshape(featureMap, [h, w]);
            featureMap = featureMap.';

        end
        
        
        % precompute height and width of feature map
        function [h, w] = get_size_featureMap(input, W, padding, stride, dilation)
            % @input: is input 
            % @W: is a weight matrix (filter)
            % @padding: zero-padding size
            % @stride: step size for traversing input
            % @dilation: factor for dilated convolution     
            % @[h, w]: height and width of feature map
                    
            % author: Dung Tran
            % date: 6/14/2019
            
            I = Conv2DLayer.get_zero_padding_input(input, padding);
            n = size(I); % n(1) is height and n(2) is width of input
            m = size(W); % m(1) is height and m(2) is width of the filter
            
            % I, W is 2D matrix
            % I is assumed to be the input after zero padding
            % output size: 
            % (InputSize - (FilterSize - 1)*Dilation + 1)/Stride
            
            h = floor((n(1) - (m(1) - 1) * dilation(1) - 1) / stride(1) + 1);  % height of feature map
            w = floor((n(2) - (m(2) - 1) * dilation(2) - 1) / stride(2) + 1);  % width of feature map

            
        end
        

        
    end
       
    % reachability analysis using star set
    
    methods
        % reachability analysis method using Stars
        % a star represent a set of images (2D matrix of h x w)
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
            c = vl_nnconv(double(input.V(:,:,:,1)), double(obj.Weights), double(obj.Bias), 'Stride', obj.Stride, 'Pad', obj.PaddingSize, 'Dilate', obj.DilationFactor);
            V = vl_nnconv(double(input.V(:,:,:,2:input.numPred + 1)), double(obj.Weights), [], 'Stride', obj.Stride, 'Pad', obj.PaddingSize, 'Dilate', obj.DilationFactor);         
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
            c = vl_nnconv(double(input.V(:,:,:,1)), double(obj.Weights), double(obj.Bias), 'Stride', obj.Stride, 'Pad', obj.PaddingSize, 'Dilate', obj.DilationFactor);
            V = vl_nnconv(double(input.V(:,:,:,2:input.numPreds + 1)), double(obj.Weights), [], 'Stride', obj.Stride, 'Pad', obj.PaddingSize, 'Dilate', obj.DilationFactor);         
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
                    error('Invalid number of input arguments, should be 1, 2, 3, 4, 5 or 6');
            end
         
            if strcmp(method, 'approx-star') || strcmp(method, 'exact-star') || strcmp(method, 'abs-dom')|| contains(method, "relax-star")
                images = obj.reach_star_multipleInputs(in_images, option);
            elseif strcmp(method, 'approx-zono')
                images = obj.reach_zono_multipleInputs(in_images, option);                
            else
                error("Unknown reachability method");
            end
            
        end
        
        
    end
    
end

