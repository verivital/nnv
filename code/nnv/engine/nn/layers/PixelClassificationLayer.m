classdef PixelClassificationLayer < handle
    % Pixel Classification Layer object to verify segmentation networks
    % Author: Dung Tran
    % Date: 4/12/2020
    
    properties
        Name = 'PixelClassificationLayer';
        Classes = [];
        OutputSize = [];
        
        NumInputs = 1;
        InputNames = {'in'};
        
    end
    
    methods
        
        function obj = PixelClassificationLayer(varargin)
            % @name: name of the layer
            % @classes: array of class
            % @outputSize: outputSize
            
            % author: Dung Tran
            % date:4/12/2020
            
            switch nargin
                
                case 3
                    
                    name = varargin{1};
                    classes = varargin{2};
                    outputSize = varargin{3};
                    
                    if ~ischar(name)
                        error('Invalid name, should be a charracter array');
                    end

                    if ~ismatrix(classes)
                        error('Invalid classes, should be a matrix');
                    end

                    if ~ismatrix(outputSize)
                        error('Invalid outputSize');
                    end           

                    if length(outputSize) ~= 3
                        error('Invalid outputSize matrix');
                    end

                    if length(classes) ~= outputSize(1, 3)
                        error('Inconsistency betwen the number of classes and the outputSize');
                    end

                    obj.Name = name;
                    obj.OutputSize = outputSize;
                    
                case 5 % used to parse the Matlab pixelClassificationLayer
                    
                    obj.Name = varargin{1};
                    classes = varargin{2};
                    obj.OutputSize = varargin{3};
                    obj.NumInputs = varargin{4};
                    obj.InputNames = varargin{5};
                    
                otherwise
                    error('Invalid number of input arguments');
            end
            
            A = categories(classes);
            B = [A; {'unknown'}; {'misclass'}]; % add two more classes for analysis
            obj.Classes = categorical(B);
            
            
            
        end
            
    end
        
        
    methods
        
        %  get an array of classes
        function classes = getClasses(obj, idxs)
            % @idxs: index array
            % author: Dung Tran
            % date: 4/22/2020
            
            n = length(idxs);
            classes = [];
            for i=1:n
                if idxs(i) > length(obj.Classes)
                    error("Invalid class index %d (input class index) > %d (maximum class index)", idxs(i), length(obj.Classes));
                end
                classes = [classes string(obj.Classes(idxs(i)))];
            end
           
        end
        
        % classified label for all pixels of an image
        function varargout = evaluate(obj,image)
            %@image: an output image before softmax layer with the size of
            %        m1 x m2 x n, where n is the number of labels needs to
            %        be classified
            % @seg_im_cat: segmentation image with class categories
            % @seg_im_id: segmentation image with class index
            
            % author: Dung Tran
            % date: 4/12/2020
            % update: 4/22/2020
                     
            n = size(image);
            if length(n)~= 3
                error('Output image should be a 3-dimensional image');
            end
            [~, y] = max(image, [], 3); 
            S = cell(n(1),n(2));
            classes = categories(obj.Classes);
            X = zeros(n(1), n(2));
            for i=1:n(1)
                for j=1:n(2)
                    S{i,j} =  classes{y(i,j)};
                    X(i,j) = y(i,j);
                end
            end
            seg_im_cat = categorical(S);
            seg_im_id = X;
            
            switch nargout
                case 2
                    varargout{1} = seg_im_id;
                    varargout{2} = seg_im_cat;
                case 1
                    varargout{1} = seg_im_id;
                otherwise
                    error("Invalid number output argument, should be 1 or 2");
            end

        end
        
        % reachability with imagestar
        function seg_im_id = reach_star_single_input(obj, IS, ~)
            % @IS: imageStar input set
            % @seg_im_id: segmentation image with class index
            
            % author: Dung Tran
            % date: 4/12/2020
            % upate: 4/22/2020, 6/1/2020
            
            h = IS.height;
            w = IS.width;
            seg_im_id = zeros(h,w);
            [im_lb, im_ub] = IS.estimateRanges; 
         
            for i=1:h
                for j=1:w
                    max_xmin = max(im_lb(i, j, :));
                    c = (im_ub(i, j, :) >= max_xmin);
                    pc = find(c);
                    if length(pc) == 1
                        seg_im_id(i,j) = pc;
                    else
                        seg_im_id(i,j) = length(obj.Classes) - 1;
                    end
                end
            end
            
            
        end
        
        % reachability with relaxed imagestar
        
        function seg_im_id = reach_relax_star_single_input(obj, IS, method, RF)
            % @IS: imageStar input set
            % @seg_im_id: segmentation image with class index
            % @method: relax-star method
            % @RF: relax factor
            
            % author: Dung Tran
            % date: 1/10/2021
            
            h = IS.height;
            w = IS.width;
            nc = IS.numChannel;
            seg_im_id = zeros(h,w);
            
            S = IS.toStar; 
            [lb, ub] = S.estimateRanges; 
            n1  = round((1-RF)*length(lb)); % number of LP need to solve
            if strcmp(method, 'relax-star-range')                
                [~,midx] = sort(ub-lb, 'descend');
                map= midx(1:n1); % neurons with optimized ranged
                xmin = S.getMins(map, 'single', 'display', 'glpk'); 
                xmax = S.getMaxs(map, 'single', 'display', 'glpk');
                lb(map) = xmin;
                ub(map) = xmax;
            elseif strcmp(method, 'relax-star-random')
                midx = randperm(length(ub), n1);
                midx = midx';             
                xmin = S.getMins(midx, 'single', 'display', 'glpk'); 
                xmax = S.getMaxs(midx, 'single', 'display', 'glpk');
                lb(midx) = xmin;
                ub(midx) = xmax;
   
            elseif strcmp(method, 'relax-star-area')
                areas = 0.5*(abs(ub).*abs(lb)); % estimated areas of triangle overapproximation at all neurons
                [~,midx] = sort(areas, 'descend');
                map= midx(1:n1); % neurons with optimized ranged
                xmin = S.getMins(map, 'single', 'display', 'glpk'); 
                xmax = S.getMaxs(map, 'single', 'display', 'glpk');
                lb(map) = xmin;
                ub(map) = xmax;
            elseif strcmp(method, 'relax-star-bound')
                N = length(ub);
                lu = [ub; abs(lb)];
                [~,midx] = sort(lu, 'descend');
                midx1 = midx(1:2*n1); % neurons with optimized ranges
                ub_idx = midx1(midx1 <= N); % neurons having upperbound optimized
                lb_idx = midx1(midx1 > N) - N;  % neurons having lowerbound optimized
                xmin = S.getMins(ub_idx, 'single', 'display', 'glpk'); 
                xmax = S.getMaxs(lb_idx, 'single', 'display', 'glpk');
                lb(lb_idx) = xmin;
                ub(ub_idx) = xmax;
            else
                error('Unknown relax-star method');
            end
                        
            im_lb = reshape(lb, [h, w, nc]);
            im_ub = reshape(ub, [h, w, nc]);

            for i=1:h
                for j=1:w
                    max_xmin = max(im_lb(i, j, :));
                    c = (im_ub(i, j, :) >= max_xmin);
                    pc = find(c);
                    if length(pc) == 1
                        seg_im_id(i,j) = pc;
                    else
                        seg_im_id(i,j) = length(obj.Classes) - 1;
                    end
                end
            end
            
            
        end
        
        % reachability with imagestar
        function seg_ims_ids = reach_star_multipleInputs(obj, in_images, option)
            % @in_images: an array of imageStar input set
            % @seg_im_id: segmentation image with class index
            % @seg_im_cat: segmentation image with class name (a
            % categorical object)
            % @mis_px_id: example of misclassified pixel class index
            
            % author: Dung Tran
            % date: 4/12/2020
            % upate: 4/22/2020
            
            
            n = length(in_images);
            seg_ims_ids = cell(n, 1);
            if strcmp(option, 'parallel')
                parfor i=1:n
                    seg_ims_ids{i} = obj.reach_star_single_input(in_images(i));
                end
            elseif strcmp(option, 'single') || isempty(option)
                for i=1:n
                    seg_ims_ids{i} = obj.reach_star_single_input(in_images(i));
                end
            else
                error('Unknown computation option');

            end           
            
            
        end
        
        
        % reachability with imagestar
        function seg_ims_ids = reach_relax_star_multipleInputs(obj, in_images, method, RF, option)
            % @in_images: an array of imageStar input set
            % @method: relax-star method 'relax-star-area',
            % 'relax-star-random', 'relax-star-bound', 'relax-star-range'
            % @seg_im_id: segmentation image with class index
                       
            % author: Dung Tran
            % date: 1/10/2020
            % upate: 1/10/2020
            
            
            n = length(in_images);
            seg_ims_ids = cell(n, 1);
            if strcmp(option, 'parallel')
                parfor i=1:n
                    seg_ims_ids{i} = obj.reach_relax_star_single_input(in_images(i), method, RF);
                end
            elseif strcmp(option, 'single') || isempty(option)
                for i=1:n
                    seg_ims_ids{i} = obj.reach_relax_star_single_input(in_images(i), method, RF);
                end
            else
                error('Unknown computation option');

            end           
            
            
        end
        
        
        % main reach method
        function seg_ims_ids = reach(varargin)
            % @in_images: an array of imageStar input set
            % @seg_im_id: segmentation image with class index
            
            % author: Dung Tran
            % date: 4/22/2020
            % upate: 6/1/2020
            
           switch nargin
                case 5
                    obj = varargin{1};
                    in_images = varargin{2};
                    method = varargin{3};
                    option = varargin{4};
                    relaxFactor = varargin{5};
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
                    error('Invalid number of input arguments, should be 1, 2, 3 or 4');
            end
         
            if strcmp(method, 'approx-star') || strcmp(method, 'exact-star') || strcmp(method, 'abs-dom') || contains(method, 'relax-star')
                seg_ims_ids = obj.reach_star_multipleInputs(in_images, option);
            elseif strcmp(method, 'approx-zono') 
                error('NNV have not support approx-zono method yet');
            else
                error('Unknown reachability method');
            end
            
            
            
        end
     
        
    end
    
    
    methods(Static)
        % parsing method
        
        function L = parse(pixel_classification_layer)
            % @pixel_classification_layer: 
            % @L: constructed layer
                        
            % author: Dung Tran
            % date: 4/12/2020
            
            
            if ~isa(pixel_classification_layer, 'nnet.cnn.layer.PixelClassificationLayer')
                error('Input is not a Matlab nnet.cnn.layer.PixelClassificationLayer');
            end
            
            L = PixelClassificationLayer(pixel_classification_layer.Name, pixel_classification_layer.Classes, pixel_classification_layer.OutputSize, pixel_classification_layer.NumInputs, pixel_classification_layer.InputNames);
            
            fprintf('\nParsing a Matlab pixel classification layer is done successfully');
            
        end
        
        
    end
end

