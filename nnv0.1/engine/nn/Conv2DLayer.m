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
        % Hyperparameters
        FilterSize = []; % height and width of filters
        NumChannels = 'auto';
        NumFilters = 0; % number of filters
        Stride = [1 1]; % step size for traversing input
        DilationFactor = [1 1]; % factor for dilated convolution
        PaddingMode = 'manual';
        PaddingSize = [0 0 0 0]; % size of padding
        
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
            
            obj.NumFilters = NumFilters;
             
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
            
            obj.Weights = weights;
            obj.Bias = biases;
            
        end
        
    end
end

