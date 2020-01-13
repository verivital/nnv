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
            % @input: high-dimensional array, for example, input(:, :, :), 
            % @y: output
            
            % author: Dung Tran
            % date: 6/17/2019
            
            % @y: high-dimensional array (output volume), depth of output = number of filters
            
            % author: Dung Tran
            % date: 12/10/2018
            % update: 7/26/2019
           
            y = vl_nnpool(double(input), obj.PoolSize, 'Stride', obj.Stride, 'Pad', obj.PaddingSize, 'Method', 'max');         
                   
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
            new_V(:,:,input.numChannel, input.numPred+1) = zeros(h,w);
            
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
                            new_V(:,:,k,p) = input.V(max_ind(1), max_ind(2), k, p);
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
        
        function images = reach_star_exact(obj, in_image)
            % @in_image: an ImageStar input set
            % @option: = 'single' single core for computation
            %          = 'parallel' multiple cores for parallel computation
            % @images: an set of imagestar with number of channels = obj.NumFilters
            
            % author: Dung Tran
            % date: 6/17/2019
            % updates: 7/25/2019
            
            if ~isa(in_image, 'ImageStar')
                error('The input is not an ImageStar object');
            end
            
            startPoints = obj.get_startPoints(in_image.V(:,:,1,1)); % get starpoints
            [h, w] = obj.get_size_maxMap(in_image.V(:,:,1,1)); % size of maxMap           
            image = in_image;   
            % check max_id first           
            max_index = cell(h, w, in_image.numChannel);
            maxMap_basis_V(:,:,in_image.numChannel, in_image.numPred+1) = zeros(h,w); % pre-allocate basis matrix for the maxmap
            split_pos = [];
            
            % compute max_index and split_index when applying maxpooling operation
            for k=1:in_image.numChannel
                for i=1:h
                    for j=1:w
                        max_index{i, j, k}  = in_image.get_localMax_index(startPoints{i,j},obj.PoolSize, k);                       
                        % construct the basis image for the maxMap
                        if size(max_index{i, j, k}, 1) == 1
                            maxMap_basis_V(i,j,k, :) = in_image.V(max_index{i, j, k}(1), max_index{i, j, k}(2), k, :);
                        else
                            split_pos = [split_pos; [i j k]];
                        end
                    end
                end
            end
            
            
            n = size(split_pos, 1);
            fprintf('\nThere are splits happened at %d local regions when computing the exact max maps', n);
            images = ImageStar(maxMap_basis_V, in_image.C, in_image.d, in_image.pred_lb, in_image.pred_ub);
            if n > 0         
                for i=1:n
                    m1 = length(images);           
                    images = ImageStar.stepSplitMultipleInputs(images, in_image, split_pos(i, :), max_index{split_pos(i, 1), split_pos(i, 2), split_pos(i, 3)}, []);
                    m2 = length(images);
                    fprintf('\nSplit %d images into %d images', m1, m2);
                end
            end
                                          
        end
        
        % reach exact star multiple inputs
        function IS = reach_star_exact_multipleInputs(obj, in_images, option)
            % @in_images: an array of imagestar input sets
            % images: an array of imagestar output sets
            % option: '[]' or 'parallel'
            
            % author: Dung Tran
            % date: 7/24/2019
            
            
            n = length(in_images);
            IS = [];
            if strcmp(option, 'parallel')
                
                parfor i=1:n
                    IS = [IS obj.reach_star_exact(in_images(i))];
                end
                
            elseif isempty(option) || strcmp(option, 'single')
                for i=1:n
                    IS = [IS obj.reach_star_exact(in_images(i))];
                end
            else
                error('Unknown computation option');
            end
            
            
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
                        
            % compute max_index when applying maxpooling operation
            for k=1:in_image.numChannel
                for i=1:h
                    for j=1:w
                        max_index{i, j, k}  = in_image.get_localMax_index(startPoints{i,j},obj.PoolSize, k);     
                    end
                end
            end
            
            % construct an over-approximate imagestar reachable set
            
            % compute new number of predicate
            np = in_image.numPred;
            l = 0;
            for k=1:in_image.numChannel
                for i=1:h
                    for j=1:w
                        max_id = max_index{i,j,k};
                        if size(max_id,1) > 1
                            np = np + 1;
                            l  = l + 1;
                        end                       
                    end
                end
            end
            
            fprintf('\n%d new variables are introduced', l);
                                   
            % update new basis matrix
            new_V(:,:,in_image.numChannel,np+1) = zeros(h,w);
            new_pred_index = 0;
            for k=1:in_image.numChannel
                for i=1:h
                    for j=1:w
                        max_id = max_index{i,j,k};
                        if size(max_id, 1) == 1                            
                            for p=1:in_image.numPred + 1
                                new_V(i,j,k, p) = in_image.V(max_id(1),max_id(2),k, p);
                            end
                        else
                            % adding new predicate variable
                            new_V(i,j,k,1) = 0;
                            new_pred_index = new_pred_index + 1;
                            new_V(i,j,k,in_image.numPred+1+new_pred_index) = 1;                            
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
                        if size(max_id,1) > 1
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
                                C2(g, 1:in_image.numPred) = in_image.V(point(1), point(2),k, 2:in_image.numPred+1);
                                C2(g, in_image.numPred + new_pred_index) = -1;
                                d2(g) = -in_image.V(point(1),point(2),k,1);                                
                            end
                            
                            C = [C1; C2];
                            d = [d1; d2];
                            
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
        
        % reach approx-star with multiple inputs
        function IS = reach_star_approx_multipleInputs(obj, in_images, option)
            % @in_images: an array of imagestar input sets
            % images: an array of imagestar output sets
            % option: '[]' or 'parallel'
            
            % author: Dung Tran
            % date: 7/24/2019
            
            
            n = length(in_images);
            IS(n) = ImageStar;
            if strcmp(option, 'parallel')
                
                parfor i=1:n
                    IS(i) = obj.reach_star_approx(in_images(i));
                end
                
            elseif isempty(option) || strcmp(option, 'single')
                for i=1:n
                    IS(i) = obj.reach_star_approx(in_images(i));
                end
            else
                error('Unknown computation option');
            end
            
            
        end
        
        
        % general functio for reachability analysis
        function IS = reach(varargin)
            % @in_image: an input imagestar
            % @image: output set
            % @option: = 'single' or 'parallel' 
            
            % author: Dung Tran
            % date: 6/26/2019
             
            switch nargin
                
                case 4
                    obj = varargin{1};
                    in_images = varargin{2};
                    method = varargin{3};
                    option = varargin{4};
                
                case 3
                    obj = varargin{1};
                    in_images = varargin{2};
                    method = varargin{3};
                    option = [];
                
                otherwise
                    error('Invalid number of input arguments (should be 2 or 3)');
            end
            
            if strcmp(method, 'approx-star')
                IS = obj.reach_star_approx_multipleInputs(in_images, option);
            elseif strcmp(method, 'exact-star')
                IS = obj.reach_star_exact_multipleInputs(in_images, option);
            elseif strcmp(method, 'abs-dom')
                % abs-domain works similarly to approx-star method for max
                % pooling layer
                IS = obj.reach_star_approx_multipleInputs(in_images, option);
            elseif strcmp(method, 'approx-zono')
                IS = obj.reach_zono_multipleInputs(in_images, option);
            end
   
            
        end
        
        
    end
    
    methods
        % reachability analysis using ImageZono
        % this is an over-approximation method
        function image = reach_zono(obj, in_image)
            % @in_image: an ImageZono input
            % @image: an ImageZono output
            
            % author: Dung Tran
            % date: 1/5/2020
            
            
            if ~isa(in_image, 'ImageZono')
                error('Input is not an ImageZono');
            end
            
            lb = in_image.lb_image;
            ub = in_image.ub_image;
            
            new_lb = vl_nnpool(-lb, obj.PoolSize, 'Stride', obj.Stride, 'Pad', obj.PaddingSize, 'Method', 'max');         
            new_ub = vl_nnpool(ub, obj.PoolSize, 'Stride', obj.Stride, 'Pad', obj.PaddingSize, 'Method', 'max');
            
            image = ImageZono(-new_lb, new_ub);
            
        end
        
        % handle multiple inputs
        function images = reach_zono_multipleInputs(obj, in_images, option)
            % @in_images: an array of zonotopes
            % @option: = 'parallel' or 'single'
            % @images: output set
            
            % author: Dung Tran
            % date: 1/6/2020
            
            
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
        
        
    end
    
    
    methods(Static)
       
        
        
        % parse a trained averagePooling2dLayer from matlab
        function L = parse(max_Pooling_2d_Layer)
            % @average_Pooling_2d_Layer: a average pooling 2d layer from matlab deep
            % neural network tool box
            % @L : a AveragePooling2DLayer obj for reachability analysis purpose
            
            % author: Dung Tran
            % date: 6/17/2019
            
            
            if ~isa(max_Pooling_2d_Layer, 'nnet.cnn.layer.MaxPooling2DLayer')
                error('Input is not a Matlab nnet.cnn.layer.MaxPooling2DLayer class');
            end
            
            L = MaxPooling2DLayer(max_Pooling_2d_Layer.Name, max_Pooling_2d_Layer.PoolSize, max_Pooling_2d_Layer.Stride, max_Pooling_2d_Layer.PaddingSize);
            fprintf('\nParsing a Matlab max pooling 2d layer is done successfully');
            
        end
        
    end
    
    
    
end

