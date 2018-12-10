classdef Conv2DLayer < handle
    % The convolutional 2D layer class
    %   Contain constructor and reachability analysis methods
    % Main references:
    % 1) An intuitive explanation of convolutional neural networks: 
    %    https://ujjwalkarn.me/2016/08/11/intuitive-explanation-convnets/
    % 2) More detail about mathematical background of CNN
    %    http://cs231n.github.io/convolutional-networks/
    % 3) Matlab implementation of Convolution2DLayer (for training and evaluating purpose)
    %    https://www.mathworks.com/help/deeplearning/ug/layers-of-a-convolutional-neural-network.html#ftn.d120e2890
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
    end
    
    
    % setting hyperparameters method
    methods
        
        % constructor of the class
        function obj = Conv2DLayer(filterSize, numFilters)
            % @filterSize: height and width of filters
            % @numFilters: number of filters
            
            % author: Dung Tran
            % date: 12/5/2018    
            % update: 
            
            [n, m] = size(filterSize);
            
            if n ~= 1
                error('Filter size should has one row');
            end
            
            if m == 1
                obj.FilterSize = [filterSize filterSize];
            elseif m == 2
                obj.FilterSize = [filterSize(1) filterSize(2)];
            elseif m ~= 1 && m ~=2
                error('Invalid filter size');
            end
            
            if numFilters < 1
                error('Number of Filters is at least one');
            end
            
            obj.NumFilters = numFilters;
             
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
            
            if length(w) ~= 4
                error('Invalid weights array');
            end
            if length(b) ~= 3
                error('Invalid biases array');
            end
            
            if w(4) ~= b(3)
                error('Inconsistency between weights and biases');
            end
            
            if w(4) ~= obj.NumFilters
                fprintf('\nInconsistency between weights and biases and number of filters parameter');
                fprintf('\nNumber of filters is updated from %d to %d', obj.NumFilters, w(4));
            end
            obj.NumFilters = w(4);
            
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
    
    
    methods(Static)
        
        % parse a trained convolutional2dLayer from matlab
        function L = parse(conv2dLayer)
            % @conv2dLayer: a convolutional 2d layer from matlab deep
            % neural network tool box
            % @L : a Cov2DLayer for reachability analysis purpose
            
            % author: Dung Tran
            % date: 12/5/2018
            
            
            if ~isa(cov2dLayer, 'convolution2dLayer')
                error('Input is not a Matlab convolution2dLayer class');
            end
            
            L = Conv2DLayer(conv2dLayer.FilterSize, conv2dLayer.NumFilters);
            L.set_stride(conv2dLayer.Stride);
            L.set_dilation(conv2dLayer.DilationFactor);
            L.set_padding(conv2dLayer.PaddingSize);
            L.set_weights_biases(conv2dLayer.Weights, conv2dLayer.Bias);
            
            fprintf('Parsing a Matlab convolutional 2d layer is done successfully');
            
        end
        
        % parse input image with padding
        function I = get_input(input, paddingSize)
            % @input: an array of input image, 1 or high-dimensional array
            % @paddingSize: paddingSize to construct the new input I
            % @I: the new input array affer applying padding
            
            % author: Dung Tran
            % date: 12/10/2018
            
            n = size(input);
            m = size(paddingSize);
            if length(m)~= 2 || m(1) ~= 1 || m(2) ~= 4
                error('Invalid PaddingSize');
            else
                t = paddingSize(1);
                b = paddingSize(2);
                l = paddingSize(3);
                r = paddingSize(4);
            end
            
            if length(n) == 2 
                % input volume has only one channel                
                h = n(1); % height of input
                w = n(2); % width of input 
                      
                I = zeros(t + h + b, l + w + r);
                I(t+1:t+h,l+1:l+w) = input; % new constructed input volume
                
            elseif length(n) > 2
                % input volume may has more than one channel
                
                % input volume has only one channel                
                h = n(1); % height of input
                w = n(2); % width of input 
                d = n(3); % depth of input (number of channel)
                
                for i=1:d
                    I(:,:,i) = zeros(t + h + b, l + w + r);
                    I(t+1:t+h, l+1:l+w, i) = input(:, :, i); % new constructed input volume
                    
                end              
                
                
            else
                error('Invalid input');
            end
                
            
            
            
        end
        
    end
    
end

