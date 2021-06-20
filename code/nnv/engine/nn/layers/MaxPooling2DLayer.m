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
        NumInputs = 1; 
        InputNames = {'in'};       
        HasUnpoolingOutputs = 0;
        NumOutputs = 1; % by default
        OutputNames = {'out'}; %by default
        
        
        % Hyperparameters
        PoolSize = [2 2];
        Stride = [1 1]; % step size for traversing input
        PaddingMode = 'manual'; 
        PaddingSize = [0 0 0 0]; % size of padding [t b l r] for nonnegative integers
        
        % used for maxunpooling
        MaxIndx = [];
        InputSize = [];
        
    end
    
    
    % setting hyperparameters method
    methods
        
        % constructor of the class
        function obj = MaxPooling2DLayer(varargin)           
            % author: Dung Tran
            % date: 6/17/2019    
            % update: 
            
            switch nargin
                
                case 9
                    
                    name = varargin{1};
                    poolSize = varargin{2};
                    stride = varargin{3};
                    paddingSize = varargin{4};
                    hasUnpoolingOutputs = varargin{5};
                    numInputs = varargin{6};
                    inputNames = varargin{7};
                    numOutputs = varargin{8};
                    outputNames = varargin{9};
                    
                    
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
                    
                    if hasUnpoolingOutputs < 0
                        error('Invalid HasUnpoolingOutputs parameter');
                    else
                        obj.HasUnpoolingOutputs = hasUnpoolingOutputs;
                    end
                    
                    if numInputs < 1
                        error('Invalid number of inputs');
                    else
                        obj.NumInputs = numInputs;
                    end
                    
                    if ~iscell(inputNames)
                        error('InputNames should be a cell');
                    else
                        obj.InputNames = inputNames;
                    end
                    
                    if numOutputs < 1
                        error('Invalid number of outputs');
                    else
                        obj.NumOutputs = numOutputs;
                    end
                    
                    if ~iscell(outputNames)
                        error('Invalid output Names, should be a cell');
                    else
                        obj.OutputNames = outputNames;
                    end
                    
                
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
                    error('Invalid number of inputs (should be 0, 3, 4, or 9)');
                                 
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
        
        function y = evaluate(varargin)
            % @input: high-dimensional array, for example, input(:, :, :), 
            % @option = 'cnn' or 'segnet'
            % @y: output
            
            % author: Dung Tran
            % date: 6/17/2019
            % update: 4/19/2020: evaluate method for segmentation network
            
            % @y: high-dimensional array (output volume), depth of output = number of filters
            
            % author: Dung Tran
            % date: 12/10/2018
            % update: 7/26/2019
            
            switch nargin
                case 2
                    obj = varargin{1};
                    input = varargin{2};
                    option = 'cnn';
                case 3
                    obj = varargin{1};
                    input = varargin{2};
                    option = varargin{3};
                otherwise
                    error('Invalid number of input');
            end
            
            if strcmp(option, 'cnn')
                
                y = vl_nnpool(double(input), obj.PoolSize, 'Stride', obj.Stride, 'Pad', obj.PaddingSize, 'Method', 'max');
                
            elseif strcmp(option, 'segnet')
                % require Matlab version 2019b
                
                dlX = dlarray(input, 'SSC'); % convert the numeric array to dlarray
                [dlY, indx, inputSize] = maxpool(dlX, obj.PoolSize, 'Stride', obj.Stride, 'Pad', obj.PaddingSize);
                y = extractdata(dlY);
                obj.MaxIndx = indx;
                obj.InputSize = inputSize;

            
            else
                error('Unknown option for max pooling evaluation');
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
        
        % get zero-padding image star
        function pad_ims = get_zero_padding_imageStar(obj, ims)
            % @ims: an image star input set
            % @pad_ims: an zero-padding image star set
            
            % author: Dung Tran
            % date: 4/23/2020
            
            if sum(obj.PaddingSize) == 0
                pad_ims = ims;
            else
                c = obj.get_zero_padding_input(ims.V(:,:,:,1));
                k = size(c);
                n = ims.numPred;
                V1(:,:,:,n+1) = zeros(k);
                for i=1:n
                    V1(:,:,:,i+1) = obj.get_zero_padding_input(ims.V(:,:,:,i+1));
                end
                V1(:,:,:,1) = c;
                if ~isempty(ims.im_lb)
                    new_im_lb = obj.get_zero_padding_input(ims.im_lb);
                    new_im_ub = obj.get_zero_padding_input(ims.im_ub);
                else
                    new_im_lb = [];
                    new_im_ub = [];
                end
                pad_ims = ImageStar(V1, ims.C, ims.d, ims.pred_lb, ims.pred_ub, new_im_lb, new_im_ub);
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
 
        function images = reach_star_exact(varargin)
            % @in_image: an ImageStar input set
            % @option: = 'single' single core for computation
            %          = 'parallel' multiple cores for parallel computation
            % @images: an set of imagestar with number of channels = obj.NumFilters
            % @dis_opt: display option: [] or 'display'
            
            % author: Dung Tran
            % date: 6/17/2019
            % updates: 7/25/2019, 4/28/2020
            
            % update: 7/15/2020: add display option
            
            switch nargin
                case 2
                    obj = varargin{1};
                    in_image = varargin{2};
                    dis_opt = [];
                    lp_solver = 'linprog';
                case 3
                    obj = varargin{1};
                    in_image = varargin{2};
                    dis_opt = varargin{3};
                    lp_solver = 'linprog';
                case 4
                    obj = varargin{1};
                    in_image = varargin{2};
                    dis_opt = varargin{3};
                    lp_solver = varargin{4};
                otherwise
                    error('Invalid number of input arguments, should be 1, 2 or 3');       
            end
            
            if ~isa(in_image, 'ImageStar')
                error('The input is not an ImageStar object');
            end
            
            startPoints = obj.get_startPoints(in_image.V(:,:,1,1)); % get starpoints
            [h, w] = obj.get_size_maxMap(in_image.V(:,:,1,1)); % size of maxMap           
            
            % padding in_image star
            pad_image = obj.get_zero_padding_imageStar(in_image);
            
            % check max_id first           
            max_index = cell(h, w, pad_image.numChannel);
            maxMap_basis_V(:,:,pad_image.numChannel, pad_image.numPred+1) = zeros(h,w); % pre-allocate basis matrix for the maxmap
            split_pos = [];
            
            % compute max_index and split_index when applying maxpooling operation
            maxidx = cell(h, w, pad_image.numChannel);
            for k=1:pad_image.numChannel
                for i=1:h
                    for j=1:w
                        max_index{i, j, k}  = pad_image.get_localMax_index(startPoints{i,j},obj.PoolSize, k, lp_solver);                       
                        % construct the basis image for the maxMap
                        if size(max_index{i, j, k}, 1) == 1
                            maxMap_basis_V(i,j,k, :) = pad_image.V(max_index{i, j, k}(1), max_index{i, j, k}(2), k, :);
                            maxidx{i, j, k} = max_index{i, j, k};
                        else
                            split_pos = [split_pos; [i j k]];
                        end
                    end
                end
            end
                        
            n = size(split_pos, 1);
            if strcmp(dis_opt, 'display')
                fprintf('\nThere are splits happened at %d local regions when computing the exact max maps', n);
            end
            images = ImageStar(maxMap_basis_V, pad_image.C, pad_image.d, pad_image.pred_lb, pad_image.pred_ub);
            images.addMaxIdx(obj.Name, maxidx);
            images.addInputSize(obj.Name, [pad_image.height pad_image.width]);
            if n > 0         
                for i=1:n
                    m1 = length(images);           
                    images = obj.stepSplitMultipleInputs(images, pad_image, split_pos(i, :, :), max_index{split_pos(i, 1), split_pos(i, 2), split_pos(i, 3)}, []);
                    m2 = length(images);
                    if strcmp(dis_opt, 'display')
                        fprintf('\nSplit %d images into %d images', m1, m2);
                    end
                end
            end
                        
                                          
        end
        
        % step split of an image star
        % a single in_image can be splitted into several images in the
        % exact max pooling operation
        function images = stepSplit(varargin)
            % @in_image: the current maxMap ImageStar
            % @ori_image: the original ImageStar to compute the maxMap 
            % @pos: local position of the maxMap where splits may occur
            % @split_index: indexes of local pixels where splits occur
            
            % author: Dung Tran
            % date: 7/25/2019
            % update: 7/16/2020: add lp_solver option
            
            switch nargin
                case 5
                    obj = varargin{1};
                    in_image = varargin{2};
                    ori_image = varargin{3};
                    pos = varargin{4};
                    split_index = varargin{5};
                    lp_solver = 'linprog';
                case 6 
                    obj = varargin{1};
                    in_image = varargin{2};
                    ori_image = varargin{3};
                    pos = varargin{4};
                    split_index = varargin{5};
                    lp_solver = varargin{6};
                otherwise
                    error('Invalid number of input arguments, should be 4 or 5');
            end
            
            
            if ~isa(in_image, 'ImageStar')
                error('input maxMap is not an ImageStar');
            end
            if ~isa(ori_image, 'ImageStar')
                error('reference image is not an ImageStar');
            end
            
            n = size(split_index);
            if n(2) ~= 3 || n(1) < 1
                error('Invalid split index, it should have 3 columns and at least 1 row');
            end
            
                        
            images = [];
            for i=1:n(1)               
                center = split_index(i, :, :);
                others = split_index;
                others(i,:,:) = [];     
                [new_C, new_d] = ImageStar.isMax(in_image, ori_image, center, others, lp_solver);                
                if ~isempty(new_C) && ~isempty(new_d)                    
                    V = in_image.V;
                    V(pos(1), pos(2), pos(3), :) = ori_image.V(center(1), center(2), center(3), :);
                    im = ImageStar(V, new_C, new_d, in_image.pred_lb, in_image.pred_ub, in_image.im_lb, in_image.im_ub);
                    im.MaxIdxs = in_image.MaxIdxs; % inherit the max indexes from the previous intermediate imagestar
                    im.InputSizes = in_image.InputSizes; % inherit the InputSize from the preivous intermediate imageStar
                    im.updateMaxIdx(obj.Name, center, pos);
                    images = [images im];
                end
            end
           
        end
        
        
        % step split for multiple image stars
        % a single in_image can be splitted into several images in the
        % exact max pooling operation
        function images = stepSplitMultipleInputs(varargin)
            % @in_image: the current maxMap ImageStar
            % @ori_image: the original ImageStar to compute the maxMap 
            % @pos: local position of the maxMap where splits may occur
            % @split_index: indexes of local pixels where splits occur
            % @option: = [] or 'parallel'
            
            % author: Dung Tran
            % date: 7/25/2019
            % update: 7/15/2020: add lp_solver option
            
            switch nargin
                case 6
                    obj = varargin{1};
                    in_images = varargin{2};
                    ori_image = varargin{3};
                    pos = varargin{4};
                    split_index = varargin{5};
                    option = varargin{6};
                    lp_solver = 'linprog';
                case 7
                    obj = varargin{1};
                    in_images = varargin{2};
                    ori_image = varargin{3};
                    pos = varargin{4};
                    split_index = varargin{5};
                    option = varargin{6};
                    lp_solver = varargin{7};
                otherwise
                    error('Invalid number of input arguments, shoule be 5 or 6');
            end
            
            
            n = length(in_images);
            images = [];
            if strcmp(option, 'parallel')
                parfor i=1:n
                    images = [images obj.stepSplit(in_images(i), ori_image, pos, split_index, lp_solver)];
                end
            elseif isempty(option) || strcmp(option, 'single')
                for i=1:n
                    images = [images obj.stepSplit(in_images(i), ori_image, pos, split_index, lp_solver)];
                end
            else 
                error('Unknown computation option');
            end       
            
        end
        
        
        % reach exact star multiple inputs
        function IS = reach_star_exact_multipleInputs(varargin)
            % @in_images: an array of imagestar input sets
            % images: an array of imagestar output sets
            % option: '[]' or 'parallel'
            % dis_opt: display option = [] or 'display'
            
            % author: Dung Tran
            % date: 7/24/2019
            % update: 7/15/2020: add display option
            %         7/16/2020: add lp_solver option
            
            switch nargin
                case 3
                    obj = varargin{1};
                    in_images = varargin{2};
                    option = varargin{3};
                    dis_opt = [];
                    lp_solver = 'linprog';
                case 4
                    obj = varargin{1};
                    in_images = varargin{2};
                    option = varargin{3};
                    dis_opt = varargin{4};
                    lp_solver = 'linprog';
                case 5
                    obj = varargin{1};
                    in_images = varargin{2};
                    option = varargin{3};
                    dis_opt = varargin{4};
                    lp_solver = varargin{5};
                otherwise
                    error('Invalid number of input arguments, should be 2, 3 or 4');       
            end
            
            n = length(in_images);
            IS = [];
            if strcmp(option, 'parallel')
                
                parfor i=1:n
                    IS = [IS obj.reach_star_exact(in_images(i), dis_opt, lp_solver)];
                end
                
            elseif isempty(option) || strcmp(option, 'single')
                for i=1:n
                    IS = [IS obj.reach_star_exact(in_images(i), dis_opt, lp_solver)];
                end
            else
                error('Unknown computation option');
            end
                
        end
        
               
        
        
        
    end
    
    
    methods % Over approximate reachability analysis using imagestar
                
        % reach star approx
        function image = reach_star_approx(varargin)
            % @in_image: input imageStar set
            % @channel_id: channel index
            % @image: output imagestar, an over-approximation of the exact
            % output set
            
            % author: Dung Tran
            % date: 6/24/2019
            % update: 7/15/2020: add display option
            %         7/16/2020: add lp_solver option
            
            switch nargin
                case 2
                    obj = varargin{1};
                    in_image = varargin{2};
                    dis_opt = [];
                    lp_solver = 'linprog';
                case 3
                    obj = varargin{1};
                    in_image = varargin{2};
                    dis_opt = varargin{3};
                    lp_solver = 'linprog';
                case 4
                    obj = varargin{1};
                    in_image = varargin{2};
                    dis_opt = varargin{3};
                    lp_solver = varargin{4};
                otherwise
                    error('Invalid number of input arguments, should be 1, 2, or 3');       
            end
            
            if ~isa(in_image, 'ImageStar')
                error('Input image is not an imagestar');
            end
            
            [h, w] = obj.get_size_maxMap(in_image.V(:,:,1,1));
            startPoints = obj.get_startPoints(in_image.V(:,:,1,1));
            max_index = cell(h, w, in_image.numChannel);
            
            % padding in_image star
            pad_image = obj.get_zero_padding_imageStar(in_image);
                        
            % compute max_index when applying maxpooling operation
            % compute new number of predicate
            np = pad_image.numPred;
            l = 0;
            for k=1:pad_image.numChannel
                for i=1:h
                    for j=1:w
                        max_index{i, j, k}  = pad_image.get_localMax_index(startPoints{i,j},obj.PoolSize, k, lp_solver);
                        max_id = max_index{i,j,k};
                        if size(max_id,1) > 1
                            np = np + 1;
                            l  = l + 1;
                        end   
                    end
                end
            end
            
            % construct an over-approximate imagestar reachable set
            if strcmp(dis_opt, 'display')
                fprintf('\n%d new variables are introduced\n', l);
            end                   
            % update new basis matrix
            new_V(:,:,pad_image.numChannel,np+1) = zeros(h,w);
            new_pred_index = 0;
            for k=1:pad_image.numChannel
                for i=1:h
                    for j=1:w
                        max_id = max_index{i,j,k};
                        if size(max_id, 1) == 1                            
                            for p=1:pad_image.numPred + 1
                                new_V(i,j,k, p) = pad_image.V(max_id(1),max_id(2),k, p);
                            end
                        else
                            % adding new predicate variable
                            new_V(i,j,k,1) = 0;
                            new_pred_index = new_pred_index + 1;
                            new_V(i,j,k,pad_image.numPred+1+new_pred_index) = 1;                            
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
            for k=1:pad_image.numChannel
                for i=1:h
                    for j=1:w
                        max_id = max_index{i,j,k};
                        if size(max_id,1) > 1
                            % construct new set of constraints here
                            new_pred_index = new_pred_index + 1;                           
                            startpoint = startPoints{i,j};
                            points = pad_image.get_localPoints(startpoint, obj.PoolSize);
                            C1 = zeros(1, np);
                            C1(pad_image.numPred + new_pred_index) = 1;
                            [lb, ub] = pad_image.get_localBound(startpoint, obj.PoolSize, k, lp_solver);
                            new_pred_lb(new_pred_index) = lb;
                            new_pred_ub(new_pred_index) = ub;
                            d1 = ub;                           
                            C2 = zeros(N, np);
                            d2 = zeros(N, 1);
                            for g=1:N
                                point = points(g,:);
                                % add new predicate constraint: xi - y <= 0
                                C2(g, 1:pad_image.numPred) = pad_image.V(point(1), point(2),k, 2:pad_image.numPred+1);
                                C2(g, pad_image.numPred + new_pred_index) = -1;
                                d2(g) = -pad_image.V(point(1),point(2),k,1);                                
                            end
                            
                            C = [C1; C2];
                            d = [d1; d2];
                            
                            new_C((new_pred_index-1)*(N+1) + 1:new_pred_index*(N+1), :) = C;
                            new_d((new_pred_index-1)*(N+1) + 1:new_pred_index*(N+1)) = d;
                        end
                    end
                end
            end
            
            n = size(pad_image.C, 1);
            C = [pad_image.C zeros(n, new_pred_index)];
            new_C = [C; new_C];
            new_d = [pad_image.d; new_d];
            new_pred_lb = [pad_image.pred_lb; new_pred_lb];
            new_pred_ub = [pad_image.pred_ub; new_pred_ub];
            
            image = ImageStar(new_V, new_C, new_d, new_pred_lb, new_pred_ub);
            image.addInputSize(obj.Name, [pad_image.height pad_image.width]);
            image.addMaxIdx(obj.Name, max_index);           
        end
        
        % reach approx-star with multiple inputs
        function IS = reach_star_approx_multipleInputs(varargin)
            % @in_images: an array of imagestar input sets
            % images: an array of imagestar output sets
            % option: '[]' or 'parallel'
            % dis_opt: = [] -> no display, 'display' -> display
            
            % author: Dung Tran
            % date: 7/24/2019
            % update: 7/15/2020: add display option
            %         7/16/2020: add lp_solver option
            
            switch nargin
                
                case 3
                    obj = varargin{1};
                    in_images = varargin{2};
                    option = varargin{3};
                    dis_opt = [];
                    lp_solver = 'linprog';
                case 4
                    obj = varargin{1};
                    in_images = varargin{2};
                    option = varargin{3};
                    dis_opt = [];
                    lp_solver = 'linprog';
                case 5
                    obj = varargin{1};
                    in_images = varargin{2};
                    option = varargin{3};
                    dis_opt = varargin{4};
                    lp_solver = varargin{5};
                    
                otherwise
                    error('Invalid number of input arguments');
                    
            end
            
            n = length(in_images);
            IS(n) = ImageStar;
            if strcmp(option, 'parallel')
                
                parfor i=1:n
                    IS(i) = obj.reach_star_approx(in_images(i), dis_opt, lp_solver);
                end
                
            elseif isempty(option) || strcmp(option, 'single')
                for i=1:n
                    IS(i) = obj.reach_star_approx(in_images(i), dis_opt, lp_solver);
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
                
                 case 7
                    obj = varargin{1};
                    in_images = varargin{2};
                    method = varargin{3};
                    option = varargin{4};
                    % relaxFactor = varargin{5}; do not use
                    dis_opt = varargin{6}; 
                    lp_solver = varargin{7}; 
                
                case 6
                    obj = varargin{1};
                    in_images = varargin{2};
                    method = varargin{3};
                    option = varargin{4};
                    %relaxFactor = varargin{5}; do not use
                    dis_opt = varargin{6};
                    lp_solver = 'glpk';
                
                case 5
                    obj = varargin{1};
                    in_images = varargin{2};
                    method = varargin{3};
                    option = varargin{4};
                    %relaxFactor = varargin{5}; do not use
                    dis_opt = [];
                    lp_solver = 'glpk';
                case 4
                    obj = varargin{1};
                    in_images = varargin{2};
                    method = varargin{3};
                    option = varargin{4};
                    dis_opt = [];
                    lp_solver = 'glpk';
                
                case 3
                    obj = varargin{1};
                    in_images = varargin{2};
                    method = varargin{3};
                    option = [];
                    dis_opt = [];
                    lp_solver = 'glpk';
                
                otherwise
                    error('Invalid number of input arguments (should be 2, 3 or 4)');
            end
            
            if strcmp(method, 'approx-star') || contains(method, 'relax-star')
                IS = obj.reach_star_approx_multipleInputs(in_images, option, dis_opt, lp_solver);
            elseif strcmp(method, 'exact-star')
                IS = obj.reach_star_exact_multipleInputs(in_images, option, dis_opt, lp_solver);
            elseif strcmp(method, 'abs-dom')
                % abs-domain works similarly to approx-star method for max
                % pooling layer
                IS = obj.reach_star_approx_multipleInputs(in_images, option, dis_opt, lp_solver);
            elseif strcmp(method, 'approx-zono')
                IS = obj.reach_zono_multipleInputs(in_images, option);
            else
                error("Unknown reachability method");
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
            
            L = MaxPooling2DLayer(max_Pooling_2d_Layer.Name, max_Pooling_2d_Layer.PoolSize, max_Pooling_2d_Layer.Stride, max_Pooling_2d_Layer.PaddingSize, max_Pooling_2d_Layer.HasUnpoolingOutputs, max_Pooling_2d_Layer.NumInputs, max_Pooling_2d_Layer.InputNames, max_Pooling_2d_Layer.NumOutputs, max_Pooling_2d_Layer.OutputNames);
            fprintf('\nParsing a Matlab max pooling 2d layer is done successfully');
            
        end
        
    end
    
    
    
end

