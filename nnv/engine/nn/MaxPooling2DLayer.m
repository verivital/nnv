classdef MaxPooling2DLayer < handle
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
        Name = 'max_pooling_2d_layer';
        % Hyperparameters
        PoolSize = [2 2];
        Stride = [1 1]; % step size for traversing input
        PaddingMode = 'manual'; 
        PaddingSize = [0 0 0 0]; % size of padding [t b l r] for nonnegative integers
        
    end
    
    
    % setting hyperparameters method
    methods
        
        % constructor of the class
        function obj = MaxPooling2DLayer(varargin)           
            % author: Dung Tran
            % date: 6/17/2019    
            % update: 
            
            switch nargin
                
                case 4
                    
                    name = varargin{1};
                    poolSize = varargin{2};
                    stride = varargin{3};
                    paddingSize = varargin{4};
                    
                    if ~ischar(name)
                        error('Name is not char');
                    else
                        obj.Name = name;
                    end                    
                    
                    if size(poolSize, 1) ~= 1 || size(poolSize, 2) ~= 2
                        error('Invalid pool size');
                    else
                        obj.PoolSize = poolSize;
                    end
                    
                    if size(stride, 1) ~= 1 || size(stride, 2) ~= 2
                        error('Invalid stride');
                    else 
                        obj.Stride = stride; 
                    end
                    
                    if size(paddingSize, 1) ~= 1 || size(paddingSize, 2) ~= 4
                        error('Invalide padding size');
                    else
                        obj.PaddingSize = paddingSize;
                    end
                
                case 3
                    
                    poolSize = varargin{1};
                    stride = varargin{2};
                    paddingSize = varargin{3};
                    
                    if size(poolSize, 1) ~= 1 || size(poolSize, 2) ~= 2
                        error('Invalid pool size');
                    else
                        obj.PoolSize = poolSize;
                    end
                    
                    if size(stride, 1) ~= 1 || size(stride, 2) ~= 2
                        error('Invalid stride');
                    else 
                        obj.Stride = stride; 
                    end
                    
                    if size(paddingSize, 1) ~= 1 || size(paddingSize, 2) ~= 4
                        error('Invalide padding size');
                    else
                        obj.PaddingSize = paddingSize;
                    end
                    
                case 0
                    
                    obj.Name = 'max_pooling_2d_layer';
                    % Hyperparameters
                    obj.PoolSize = [2 2];
                    obj.Stride = [1 1]; % step size for traversing input
                    obj.PaddingMode = 'manual'; 
                    obj.PaddingSize = [0 0 0 0]; % size of padding [t b l r] for nonnegative integers
                                    
                otherwise
                    error('Invalid number of inputs (should be 0, 3 or 4)');
                                 
            end 
             
        end
        
        % set poolsize method 
        function set_poolSize(obj, poolSize)
            % @stride: stride matrix
            % author: Dung Tran
            % date: 12/5/2018
            
            [n, m] = size(poolSize);
            
            if n ~= 1
                error('poolSize matrix shoule have one row');
            end
            
            if m == 1
                obj.PoolSize = [poolSize poolSize];
            elseif m == 2
                obj.PoolSize = [poolSize(1) poolSize(2)];
            elseif m ~=1 && m ~=2 
                error('Invalid poolSize matrix');
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
                error('Invalid stride matrix');
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
                error('Invalid padding matrix');
            end
        end
        
        
    end
        
    % evaluation method
    methods
        
        function y = evaluate(obj, input)
            % @input: 2 or 3-dimensional array, for example, input(:, :, :), 
            % @y: 2 or 3-dimensional array, for example, y(:, :, :)
            
            % author: Dung Tran
            % date: 6/17/2019
            
            % @y: high-dimensional array (output volume), depth of output = number of filters
            
            % author: Dung Tran
            % date: 12/10/2018
            
            
            n = size(input);
            if length(n) == 3
                NumChannels = n(3);
            elseif length(n) == 2
                NumChannels = 1;
            else
                error('Invalid input');
            end
            
            [h, w] = obj.get_size_maxMap(input(:,:,1));   
            y(:,:,NumChannels) = zeros(h, w); % preallocate 3D matrix
            for i=1:NumChannels % number of channels
                % compute average map with i^th channels 
                y(:, :, i) = obj.compute_maxMap(input(:,:,i));
            end
            
                   
        end
        
        % parallel evaluation on an array of inputs
        function y = evaluate_parallel(obj, inputs)
            % @inputs: an array of inputs
            % @y: an array of outputs 
            
            % author: Dung Tran
            % date: 12/16/2018
            y = [];
            parfor i=1:length(inputs)
                y = [y obj.evaluate(inputs(i))];               
            end
            
            
        end
        
        
        
        % parse a trained averagePooling2dLayer from matlab
        function parse(obj, max_Pooling_2d_Layer)
            % @average_Pooling_2d_Layer: a average pooling 2d layer from matlab deep
            % neural network tool box
            % @L : a AveragePooling2DLayer obj for reachability analysis purpose
            
            % author: Dung Tran
            % date: 6/17/2019
            
            
            if ~isa(max_Pooling_2d_Layer, 'maxPooling2dLayer')
                error('Input is not a Matlab maxPooling2dLayer class');
            end
            
            obj.Name = max_Pooling_2d_Layer.Name;
            obj.PoolSize = max_Pooling_2d_Layer.PoolSize;
            obj.Stride = max_Pooling_2d_Layer.Stride;
            obj.PaddingSize = max_Pooling_2d_Layer.PaddingSize;
            fprintf('Parsing a Matlab convolutional 2d layer is done successfully');
            
        end
        
        % parse input image with padding
        function padded_I = get_zero_padding_input(obj, input)
            % @input: an array of input image
            % @paddingSize: paddingSize to construct the new input I
            % @padded_I: the new input array affer applying padding
            
            % author: Dung Tran
            % date: 6/17/2019
            
            n = size(input);        
                        
            t = obj.PaddingSize(1);
            b = obj.PaddingSize(2);
            l = obj.PaddingSize(3);
            r = obj.PaddingSize(4);
 
            
            if length(n) == 2 
                % input volume has only one channel                
                h = n(1); % height of input
                w = n(2); % width of input 
                      
                padded_I = zeros(t + h + b, l + w + r);
                padded_I(t+1:t+h,l+1:l+w) = input; % new constructed input volume
                
            elseif length(n) > 2
                % input volume may has more than one channel
                            
                h = n(1); % height of input
                w = n(2); % width of input 
                d = n(3); % depth of input (number of channel)
               
                padded_I = zeros(t + h + b, l + w + r, d); % preallocate 3D matrix
                for i=1:d
                    padded_I(t+1:t+h, l+1:l+w, i) = input(:, :, i); % new constructed input volume
                end              
                
                
            else
                error('Invalid input');
            end
                         
        end
        
        
        % compute feature map for specific input and weight
        function maxMap = compute_maxMap(obj, input)
            % @input: is input image 
            % @maxMap: maxMap of the input image
                        
            % author: Dung Tran
            % date: 6/17/2019
            
            % this maxMap computation work similar to the featureMap
            
            I = obj.get_zero_padding_input(input);
            m = obj.PoolSize; % m(1) is height and m(2) is width of the filter            
            [h, w] = obj.get_size_maxMap(input);

            % a collection of start points (the top-left corner of a square) of the region of input that is filtered
            map = obj.get_startPoints(input); 
            
            % compute feature map for each cell of map
            % do it in parallel using cpu or gpu
            % TODO: explore power of using GPU for this problem
            
            maxMap = zeros(1, h*w); % using single vector for parallel computation

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
                i0 = map{i, j}(1);
                j0 = map{i, j}(2);
                val = I(i0, j0);
                for i=i0:1:i0+m(1)-1
                    for j=j0:1:j0+m(2)-1
                        if val < I(i, j)
                            val = I(i,j);
                        end
                    end
                end
                
                maxMap(l) = val;               

            end

            maxMap = reshape(maxMap, [h, w]);
            maxMap = maxMap.';

        end
        
        
        % precompute height and width of max map
        function [h, w] = get_size_maxMap(obj, input)
            % @input: is input 
            % @[h, w]: height and width of average map
                    
            % author: Dung Tran
            % date: 6/17/2019
            
            % reference: http://cs231n.github.io/convolutional-networks/#pool
            
            I = obj.get_zero_padding_input(input);
            n = size(I); % n(1) is height and n(2) is width of input
            m = obj.PoolSize; % m(1) is height and m(2) is width of the pooling filter
            
            % I is assumed to be the input after zero padding
            % output size: 
            % (InputSize - FilterSize)/Stride + 1
            
            h = floor((n(1) - m(1) ) / obj.Stride(1) + 1);  % height of average map
            w = floor((n(2) - m(2)) / obj.Stride(2) + 1);  % width of average map

        end
        
        % get collection of start points for computing maxMap
        function startPoints = get_startPoints(obj, input)
            % @input: is input image
            % @startPoints: start points of maxMap
            
            % author: Dung Tran
            % date: 6/20/2019
            
            
            I = obj.get_zero_padding_input(input);
            m = obj.PoolSize; % m(1) is height and m(2) is width of the filter            
            [h, w] = obj.get_size_maxMap(input);

            % a collection of start points (the top-left corner of a square) of the region of input that is filtered
            startPoints = cell(h, w); 
            
            for i=1:h
                for j=1:w
                    startPoints{i, j} = zeros(1, 2);

                    if i==1
                        startPoints{i, j}(1) = 1;
                    end
                    if j==1
                        startPoints{i, j}(2) = 1;
                    end

                    if i > 1
                        startPoints{i, j}(1) = startPoints{i - 1, j}(1) + obj.Stride(1);
                    end

                    if j > 1
                        startPoints{i, j}(2) = startPoints{i, j - 1}(2) + obj.Stride(2);
                    end

                end
            end
            
        end
        
    end
        
    % reachability analysis using star set
    methods
        % reachability analysis method using Stars
        % a star represent a set of images (2D matrix of h x w)
        
        % stepMaxPooling_exact
        % compute exact imagestar of a single maxpooling step
        
        function images = stepMaxPooling_exact_singleInput(obj, image, start_point, channel_ind)
            % @image: input image, an imagestar set
            % @start_point: the top left-corner point being filter by
            % maxpooling operation
            % @channel_ind: channel index 
            
            % @images: output images after a single stepMaxPooling
            % note***: one stepMaxPooling in the worst case produces N =
            % obj.height * obj.width new images, analysis complexity is O(N)
            % number of constraints in the resulted images is O(n0 + N -
            % 1); n0 is the number of constraints of the input image
            
            % author: Dung Tran
            % date: 6/20/2019
            
            if ~isa(image, 'ImageStar')
                error('Input image is not an ImageStar');
            end
            
            
            if size(start_point, 1) ~= 2 || size(start_point, 2) ~= 1 || start_point(1) > image.height || start_point(1) < 1 || start_point(2) > image.width || start_point(2) < 1
                error('Invalid start point');
            end
            
            if channel_ind < 1 || channel_ind > image.numChannel
                error('Invalid channel index');
            end
            
            N = obj.PoolSize(1) * obj.PoolSize(2); 
            
            ind_mat = zeros(2, N);
                     
            for i=1:obj.PoolSize(1) % height of the max pool
                for j=1:obj.PoolSize(2) % width of the max pool
                    ind = (i-1)*obj.PoolSize(2) + j;
                    ind_mat(:, ind) = [start_point(1) + i-1; start_point(2) + j - 1];
                end
            end
            
            images = [];
            for i=1:N
                center = ind_mat(:,i);
                others = ind_mat;
                others(:, i) = [];
                image1 = image.isMax(center, others, channel_ind);
                images = [images image1];
            end          
            
        end
        
        % exact stepMaxPooling with multiple inputs
        function images = stepMaxPooling_exact_multipleInputs(varargin)
            % @in_images: an array of input images (an array of ImageStar)
            % @start_points: the top left-corner point being filter by
            % maxpooling operation
            % @channel_ind: channel index
            % @option: = 'parallel', use parallel computation
            %          = '[]', don't use parallel computation
            
            switch nargin
                case 4
                    obj = varargin{1};
                    in_images = varargin{2};
                    start_point = varargin{3};
                    channel_ind = varargin{4};
                    option = [];
                    
                case 5
                    obj = varargin{1};
                    in_images = varargin{2};
                    start_point = varargin{3};
                    channel_ind = varargin{4};
                    option = varargin{5};
                    if ~strcmp(option, 'parallel')
                        error('Unknow computation option');
                    end
                    
                otherwise
                    error('Invalid number of input argument, should be 4 or 5');
                
            end
            
            % do stepMaxPooling with multiple inputs
            N = length(in_images);
            images = [];
            if ~isempty(option)
                parfor i=1:N % use parallel computing
                    image1 = obj.stepMaxPooling_exact_singleInput(in_images(i), start_point, channel_ind);
                    images = [images image1];
                end
            else
                for i=1:N
                    image1 = obj.stepMaxPooling_exact_singleInput(in_images(i), start_point, channel_ind);
                    images = [images image1];
                end
            end
                    
        end
        
        
        function images = reach_star_exact(obj, input)
            % @inputs: an ImageStar input set
            % @option: = 'single' single core for computation
            %          = 'parallel' multiple cores for parallel computation
            % @S: an imagestar with number of channels = obj.NumFilters
            
            % author: Dung Tran
            % date: 6/17/2019
            
            if ~isa(input, 'ImageStar')
                error('The input is not an ImageStar object');
            end
            
            startPoints = obj.get_startPoints(input.V(:,:,1,1));
            [h, w] = size(startPoints); % this is the size of the maxMap
            
            image = input;
            
            for k=1:input.numChannel
                for i=1:h
                    for j=1:w
                        image = obj.stepMaxPooling_exact_multipleInputs(image, startPoints{i,j}', k); 
                    end
                end
            end
            
                        
            images = image;
            
            
                           
        end
    end
    
end

