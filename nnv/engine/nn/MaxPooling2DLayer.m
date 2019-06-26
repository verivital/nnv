classdef MaxPooling2DLayer < handle
    % The MaxPooling 2D layer class in CNN
    %   Contain constructor and reachability analysis methods
    % Main references:
    % 1) An intuitive explanation of convolutional neural networks: 
    %    https://ujjwalkarn.me/2016/08/11/intuitive-explanation-convnets/
    % 2) More detail about mathematical background of CNN
    %    http://cs231n.github.io/convolutional-networks/
    %    http://cs231n.github.io/convolutional-networks/#pool
    % 3) Matlab implementation of Convolution2DLayer and MaxPooling (for training and evaluating purpose)
    %    https://www.mathworks.com/help/deeplearning/ug/layers-of-a-convolutional-neural-network.html
    %    https://www.mathworks.com/help/deeplearning/ref/nnet.cnn.layer.maxpooling2dlayer.html
    
    %   Dung Tran: 6/20/2019
    
    properties
        Name = 'max_pooling_2d_layer';
        % Hyperparameters
        PoolSize = [2 2];
        Stride = [1 1]; % step size for traversing input
        PaddingMode = 'manual'; 
        PaddingSize = [0 0 0 0]; % size of padding [t b l r] for nonnegative integers
        
        NumSplits = 0; % total number of splits in the exact analysis
        
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
            fprintf('Parsing a Matlab max pooling 2d layer is done successfully');
            
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
        
        % construct an ImageStar mapMap after performing reachability analysis
        function image = construct_maxImage(obj, input)
            % @input: an input image with max_points information
            % @image: a max-ImageStar set
            
            % author: Dung Tran
            % date: 6/21/2019
            
            if ~isa(input, 'ImageStar')
                error('Input is not an ImageStar');
            end
            
            [h, w] = obj.get_size_maxMap(input.V(:,:,1,1));
            new_V(:,:,input.numPred+1,input.numChannel) = zeros(h,w);
            
            % get max points for each channel
            channel_maxPoints(:,:,input.numChannel) = zeros(3, h*w);
            for i=1:input.numChannel
                channel_maxPoints(:,:,i) = input.max_points(:, (i-1)*h*w + 1:i*h*w);
            end
            
            for p=1:input.numPred+1
                for k=1:input.numChannel
                    for i=1:h
                        for j=1:w
                            ind = (i-1)*w + j;
                            max_ind = channel_maxPoints(:,ind,k);
                            new_V(:,:,p,k) = input.V(max_ind(1), max_ind(2), p, k);
                        end
                    end
                end
            end
                        
            image = ImageStar(new_V, input.C, input.d, input.pred_lb, input.pred_ub);
            
        end
        
    end
        
    % exact reachability analysis using star set
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
            % optimize code here
            % number of feasibility problem needs to check is N
            % is there a best way to reduce this?
            for i=1:N
                center = ind_mat(:,i);
                others = ind_mat;
                others(:, i) = [];
                image1 = image.isMax(center, others, channel_ind); % optimize the isMax function
                images = [images image1];
            end
            
            obj.NumSplits = obj.NumSplits + length(images) - 1;
             
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
            obj.NumSplits = 0;
            
            for k=1:input.numChannel
                for i=1:h
                    for j=1:w
                        fprintf('\nPerforming the %d^th exact stepMaxPooling on channel %d', (i-1)*w + j, k);
                        image = obj.stepMaxPooling_exact_multipleInputs(image, startPoints{i,j}', k);
                    end
                end
            end
            
            N = length(image);
            images(1, N) = ImageStar;
            % use parallel computing for constructing the max images
            fprintf("\nConstructing the maxImage reachable sets");
            for i=1:N
                images(i) = obj.construct_maxImage(image(i));
            end
            fprintf("\nTotal number of the reachable sets after max pooling operation is %d\n", N);
                                          
        end
    end
    
    
    methods % Over approximate reachability analysis using imagestar
                
        % reach star approx
        function image = reach_star_approx(obj, in_image)
            % @in_image: input imageStar set
            % @channel_id: channel index
            % @image: output imagestar, an over-approximation of the exact
            % output set
            
            % author: Dung Tran
            % date: 6/24/2019
            
            
            if ~isa(in_image, 'ImageStar')
                error('Input image is not an imagestar');
            end
            
            [h, w] = obj.get_size_maxMap(in_image.V(:,:,1,1)); 
            startPoints = obj.get_startPoints(in_image.V(:,:,1,1));
            max_index = cell(h, w, in_image.numChannel);
            max_mat = cell(h, w, in_image.numChannel);
            
            % compute max_index when applying maxpooling operation
            for k=1:in_image.numChannel
                for i=1:h
                    for j=1:w
                        [max_index{i, j, k}, max_mat{i,j,k}] = in_image.get_localMax_index(startPoints{i,j},obj.PoolSize, k);
                    end
                end
            end
            
            % construct an over-approximate imagestar reachable set
            
            % compute new number of predicate
            np = in_image.numPred;
            for k=1:in_image.numChannel
                for i=1:h
                    for j=1:w
                        max_id = max_index{i,j,k};
                        if isempty(max_id)
                            np = np + 1;
                        end                       
                    end
                end
            end
            
            % update new basis matrix
            new_V(:,:,np+1,in_image.numChannel) = zeros(h,w);
            new_pred_index = 0;
            for k=1:in_image.numChannel
                for i=1:h
                    for j=1:w
                        max_id = max_index{i,j,k};
                        if ~isempty(max_id)                            
                            for p=1:in_image.numPred + 1
                                new_V(i,j,p,k) = in_image.V(max_id(1),max_id(2),p,k);
                            end
                        else
                            % adding new predicate variable
                            new_V(i,j,1,k) = 0;
                            new_pred_index = new_pred_index + 1;
                            new_V(i,j,in_image.numPred+1+new_pred_index,k) = 1;                            
                        end                       
                    end
                end
            end
            
            %update constraint matrix C and constraint vector d
            N = obj.PoolSize(1) * obj.PoolSize(2);            
            new_C = zeros(new_pred_index*(N + 1), np);
            new_d = zeros(new_pred_index*(N + 1), 1);
            new_pred_lb = zeros(new_pred_index, 1); % update lower bound and upper bound of new predicate variables
            new_pred_ub = zeros(new_pred_index, 1); 
            new_pred_index = 0;
            for k=1:in_image.numChannel
                for i=1:h
                    for j=1:w
                        max_id = max_index{i,j,k};
                        if isempty(max_id)
                            % construct new set of constraints here
                            new_pred_index = new_pred_index + 1;                           
                            startpoint = startPoints{i,j};
                            points = in_image.get_localPoints(startpoint, obj.PoolSize);
                            C1 = zeros(1, np);
                            C1(in_image.numPred + new_pred_index) = 1;
                            [lb, ub] = in_image.get_localBound(startpoint, obj.PoolSize, k);
                            new_pred_lb(new_pred_index) = lb;
                            new_pred_ub(new_pred_index) = ub;
                            d1 = ub;                           
                            C2 = zeros(N, np);
                            d2 = zeros(N, 1);
                            for g=1:N
                                point = points(g,:);
                                % add new predicate constraint: xi - y <= 0
                                display(points);
                                display(point);
                                C2(g, 1:in_image.numPred) = in_image.V(point(1), point(2),2:in_image.numPred+1, k);
                                C2(g, in_image.numPred + new_pred_index) = -1;
                                d2(g) = -in_image.V(point(1),point(2),1,k);                                
                            end
                            
                            C = [C1; C2];
                            d = [d1;d2];
                            
                            new_C((new_pred_index-1)*(N+1) + 1:new_pred_index*(N+1), :) = C;
                            new_d((new_pred_index-1)*(N+1) + 1:new_pred_index*(N+1)) = d;
                        end
                    end
                end
            end
            
            n = size(in_image.C, 1);
            C = [in_image.C zeros(n, new_pred_index)];
            new_C = [C; new_C];
            new_d = [in_image.d; new_d];
            new_pred_lb = [in_image.pred_lb; new_pred_lb];
            new_pred_ub = [in_image.pred_ub; new_pred_ub];
            
            image = ImageStar(new_V, new_C, new_d, new_pred_lb, new_pred_ub);
                       
        end
        
        
    end
    
end

