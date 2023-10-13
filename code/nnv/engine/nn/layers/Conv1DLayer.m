classdef Conv1DLayer < handle
    % The convolutional 1D layer class
    %   Contain constructor and reachability analysis methods
    % Main references:
    % 1) Matlab implementation of Convolution1DLayer (for training and evaluating purpose)
    %    https://www.mathworks.com/help/deeplearning/ug/layers-of-a-convolutional-neural-network.html
    %    https://www.mathworks.com/help/deeplearning/ref/nnet.cnn.layer.convolution1dlayer.html
    
    %   Neelanjana Pal: 10/25/2021
    
    properties
        Name = 'convolutional_1d_layer';
        % Hyperparameters
        FilterSize = 1; % height and width of filters
        NumChannels = 'auto';   
        NumFilters = 0; % number of filters
        Stride = 1; % step size for traversing input
        DilationFactor = 1; % factor for dilated convolution
        PaddingMode = 'manual';
        PaddingSize = [0 0]; % size of padding [l r] for nonnegative integers
        
        % Learnable Parmeters/ Used for reachability analysis
        Weights = [];
        Bias = [];
        
        %Layer parameters
        NumInputs = 1;
        InputNames = {'in'};
        NumOutputs = 1;
        OutputNames = {'out'};
    end
    
    
    % setting hyperparameters method
    methods
        
        % constructor of the class
        function obj = Conv1DLayer(varargin)                
            
            switch nargin
                
                case 10 % used for parsing a Matlab conv1dlayer 
                    
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
                    
                    if length(w) ~= 3
                        error('Invalid weights array');
                    end
                    if length(b) ~= 2
                        error('Invalid biases array');
                    end

                    if w(3) ~= b(2)
                        error('Inconsistency between filter weights and filter biases');
                    end
                    obj.NumFilters = w(3);
                    obj.NumChannels = w(2);
                    
                    obj.FilterSize = w(1);
                    
                    obj.Weights = filter_weights;
                    obj.Bias = filter_bias;
                                        
                    p = size(padding_mat);
                   
                    if  length(p) ~= 2 || p(2) ~= 2|| p(1) ~= 1
                        error('Invalid padding matrix');
                    end
                    obj.PaddingSize = padding_mat;

                    obj.Stride = stride_mat;
                                 
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
                    
                    if length(w) ~= 3
                        error('Invalid weights array');
                    end
                    if length(b) ~= 2
                        error('Invalid biases array');
                    end

                    if w(3) ~= b(2)
                        error('Inconsistency between filter weights and filter biases');
                    end
                    obj.NumFilters = w(3);
                    obj.NumChannels = w(2);
                    
                    obj.FilterSize = w(1);
                    obj.Weights = filter_weights;
                    obj.Bias = filter_bias;
                                        
                    p = size(padding_mat);
                   
                    if  length(p) ~= 2 || p(2) ~= 2|| p(1) ~= 1
                        error('Invalid padding matrix');
                    end
                    obj.PaddingSize = padding_mat;
                    
                    obj.Stride = stride_mat;
             
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
                        obj.NumChannels = w(2);
                        obj.FilterSize = w(1);
                        obj.Weights = filter_weights;
                        
                    elseif length(w) == 3
                        obj.NumFilters = w(3);
                        obj.NumChannels = w(2);
                        obj.FilterSize = w(1);
                        obj.Weights = filter_weights;
%                     elseif length(w) == 4
%                         obj.NumFilters = w(4);
%                         obj.NumChannels = w(3);
%                         obj.FilterSize = w(1);
%                         obj.Weights = filter_weights;
                    else
                        error('Invalid weight matrix');
                    end
                    
                    if length(b) ~= 2
                        error('Invalid biases array');
                    else
                        obj.Bias = filter_bias;
                    end

                    if w(3) ~= b(2)
                        error('Inconsistency between filter weights and filter biases');
                    end
                    
                    p = size(padding_mat);
                    if  length(p) ~= 2 || p(2) ~= 2|| p(1) ~= 1
                        error('Invalid padding matrix');
                    end
                    obj.PaddingSize = padding_mat;
                                               
                    obj.Stride = stride_mat;
               
                    obj.DilationFactor = dilation_mat;                   
                    
                    
                case 2
                    
                    filter_weights = varargin{1};
                    filter_bias = varargin{2};
                                        
                    obj.Name = 'convolutional_layer';                   
                    
                    w = size(filter_weights);
                    b = size(filter_bias);
                    
                    if length(w) == 2
                        obj.NumFilters = 1;
                        obj.NumChannels = w(2);
                        obj.FilterSize = w(1);
                        obj.Weights = filter_weights;
                        
                    elseif length(w) == 3
                        obj.NumFilters = w(3);
                        obj.NumChannels = w(2);
                        obj.FilterSize = w(1);
                        obj.Weights = filter_weights;
%                     elseif length(w) == 4
%                         obj.NumFilters = w(4);
%                         obj.NumChannels = w(3);
%                         obj.FilterSize = w(1);
%                         obj.Weights = filter_weights;
                    else
                        error('Invalid weight matrix');
                    end                 

                    obj.Bias = filter_bias;

                    if length(w) == 3 && w(3) ~= b(2)
                        error('Inconsistency between filter weights and filter biases');
                    end
                                        
                    % use default setting for Stride, Padding and Dilation
                    obj.Stride = 1; % step size for traversing input
                    obj.DilationFactor = 1; % factor for dilated convolution
                    obj.PaddingMode = 'manual';
                    obj.PaddingSize = [0 0]; % size of padding [l r] for nonnegative integers
        
                                    
                otherwise
                    error('Invalid number of inputs (should be 2, 5, or 6)');
                                 
            end
                    
                
            
            
            
            
             
        end

        
        % set padding 
        function set_padding(obj, padding)
            % @padding: padding matrix
            % Neelanjana Pal: 10/25/2021
            
            [n, m] = size(padding);
            
            if n ~= 1
                error('Padding matrix shoule have one row');
            end
            
            if m == 1
                obj.PaddingSize = [padding padding];
            elseif m == 2
                obj.PaddingSize = [padding(1) padding(2)];
            elseif m ~=1 && m ~=4 
                error('Invalide padding matrix');
            end
        end
        
        % set weights and biases
        function set_weights_biases(obj, weights, biases)
            % @weights: a 3-dimensional array
            %  = FilterSize-by-NumChannels-by-NumFilter
            % @biases: a 2-dimensional array
            %  = 1-by-NumFilters 
            
            w = size(weights);
            b = size(biases);
            
            if length(w) == 1
                obj.NumFilters = 1;
                obj.NumChannels = 1;
            elseif length(w) == 2
                obj.NumFilters = 1;
                obj.NumChannels = w(2);
            elseif length(w) == 3
                obj.NumFilters = w(3);
                obj.NumChannels = w(2);
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
    
    
    % evaluation method
    methods
        
        % evaluate
        function y = evaluateSequence(obj, input)
            % @input: 2-dimensional array, for example, input(:, :)
            % @y: high-dimensional array (output volume), depth of output = number of filters

            n = length(size(input));
            if n==2
                input = dlarray(input');
                y = dlconv(input, obj.Weights, obj.Bias,"Stride",obj.Stride,"DilationFactor",obj.DilationFactor,"Padding",obj.PaddingSize,"DataFormat","SC");
                y = extractdata(y)';
            end

        end
              
        
    end
    
    
    methods(Static) % parsing, get zero input pading, compute feature maps
        
        % parse a trained convolutional1dLayer from matlab
        function L = parse(conv1dLayer)
            % @conv1dLayer: a convolutional 1d layer from matlab deep
            % neural network tool box
            % @L : a Cov1DLayer for reachability analysis purpose
            
            % Neelanjana Pal: 10/25/2021
            
            if ~isa(conv1dLayer, 'nnet.cnn.layer.Convolution1DLayer')
                error('Input is not a Matlab nnet.cnn.layer.Convolution1DLayer class');
            end
            
            layer_name = conv1dLayer.Name; 
            filter_weights = conv1dLayer.Weights;
            filter_bias = conv1dLayer.Bias;
            padding_mat = conv1dLayer.PaddingSize;
            stride_mat = conv1dLayer.Stride;
            dilation_mat = conv1dLayer.DilationFactor;
            
            L = Conv1DLayer(layer_name, filter_weights, filter_bias, padding_mat, stride_mat, dilation_mat, conv1dLayer.NumInputs, conv1dLayer.InputNames, conv1dLayer.NumOutputs, conv1dLayer.OutputNames);
            
        end
        
        % parse input image with padding
        function I = get_zero_padding_input(input, paddingSize)
            % @input: an array of input image, 1 or high-dimensional array
            % @paddingSize: paddingSize to construct the new input I
            % @I: the new input array affer applying padding
            
            % author: Neelanjana Pal
            % date: 10/28/2021

            
            n = size(input);        
            if length(paddingSize) ~= 2
                error("\Invalid padding size");
            end
            
            l = paddingSize(1);
            r = paddingSize(2);
 
            
            if length(n) == 2 
                % input volume has only one channel                
                h = n(1); % height of input
                w = n(2); % width of input 
                      
                I = zeros( h, l + w + r);
                I(1:h,l+1:l+w) = input; % new constructed input volume
                
            elseif length(n) > 2
                % input volume may has more than one channel
                            
                h = n(1); % height of input
                w = n(2); % width of input 
                d = n(3); % depth of input (number of channel)
               
                I = zeros(h, l + w + r, d); % preallocate 3D matrix
                for i=1:d
                    I(1:h, l+1:l+w, i) = input(:, :, i); % new constructed input volume
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
            
            
            % author: Neelanjana Pal
            % date: 10/28/2021
            
            [h, w] = get_size_featureMap(input, W, padding, stride, dilation);

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
                    
            % author: Neelanjana Pal
            % date: 10/28/2021
            
            I = Conv1DLayer.get_zero_padding_input(input, padding);
            n = size(I); % n(1) is number of the channel, n(2) is width of input and n(3) is no of filters
            m = size(W); % m(1) is width and m(2) is number of the channel
            
            % I, W is 2D matrix
            % I is assumed to be the input after zero padding
            % output size: 
            % (InputSize - ((FilterSize - 1)*Dilation + 1)/Stride + 1)
            
            h = m(3);  % height of feature map
            w = floor((n(2) - (m(1) - 1) * dilation - 1) / stride + 1);  % width of feature map

        end
        

        
    end
       
    % reachability analysis using star set
    
    methods
        % reachability analysis method using Stars
        % a star represent a set of images (2D matrix of h x w)
        function S = reach_star_single_input(obj, input)
            % @inputs: an ImageStar input set
            % @S: an ImageStar with number of channels = obj.NumFilters
            
            if ~isa(input, 'ImageStar')
                error('The input is not an ImageStar object');
            end
            
            if input.height ~= obj.NumChannels
                error("Input set contains %d channels while the convolutional layers has %d channels", input.numChannel, obj.NumChannels);
            end

            if isempty(input.im_lb) && isempty(input.im_ub)
                layer = obj;
                c = layer.evaluateSequence(input.V(:,:,:,1));
                layer.Bias = [];
                n = size(input.V(:,:,:,2:input.numPred + 1),4);
                for i = 1:n
                    V(:,:,1,i) = layer.evaluateSequence(input.V(:,:,:,i+1));
                end
                Y = cat(4, c, V);
                S = ImageStar(Y, input.C, input.d, input.pred_lb, input.pred_ub);
            else
                lb = obj.evaluateSequence(input.im_lb);
                ub = obj.evaluateSequence(input.im_ub);
                S = ImageStar(lb,ub);
            end
                  
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
            c = vl_nnconv1d(input.V(:,:,:,1), obj.Weights, obj.Bias', obj.Stride, obj.PaddingSize, obj.DilationFactor);
            n = size(input.V(:,:,:,2:input.numPred + 1),4);
            for i = 1:n
                V(:,:,1,i) = vl_nnconv1d(input.V(:,:,:,i+1), obj.Weights, 0, obj.Stride, obj.PaddingSize, obj.DilationFactor);         
            end
            Y = cat(4, c, V);
            Z = ImageZono(Y);
            
        end
        
        
        % hangle multiple inputs
        function images = reach_star_multipleInputs(obj, in_signals, option)
            % @in_images: an array of ImageStars input set
            % @option: 
            % @images: an array of ImageStars output set
            
            n = length(in_signals);
            images(n) = ImageStar; 
            
            if strcmp(option, 'parallel')
                parfor i=1:n
                    images(i) = obj.reach_star_single_input(in_signals(i));
                end
            elseif strcmp(option, 'single') || isempty(option)
                for i=1:n
                    images(i) = obj.reach_star_single_input(in_signals(i));
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
        function images = reachSequence(varargin)
            % @inputs: an array of ImageStar or ImageZono input set
            % @S: an array of ImageStar output set
            
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

