classdef LayerNormalizationLayer < handle
    % The Layer Normalization Layer class in CNN
    %   Contain constructor and reachability analysis methods   
    %  Neelanjana Pal: 8/30/2023
    
    properties
        Name = 'LayerNormalizationLayer';
        
        NumChannels = [];
        
        % Hyperparameters 
        Epsilon = 0.00001;  % default value
        
        % Learnable parameters
        Offset = [];
        Scale = [];
        
        % layer properties
        NumInputs = 1;
        InputNames = {'in'};
        NumOutputs = 1;
        OutputNames = {'out'};
               
    end
    
    
    % setting hyperparameters method
    methods
        
        % constructor of the class
        function obj = LayerNormalizationLayer(varargin)           
            % author:Neelanjana Pal
            % date: 8/30/2023    
            % update: 
            
            
            if mod(nargin, 2) ~= 0
                error('Invalid number of arguments');
            end
            
            for i=1:nargin-1
                
                if mod(i, 2) ~= 0
                    
                    if strcmp(varargin{i}, 'Name')
                        obj.Name = varargin{i+1};
                    elseif strcmp(varargin{i}, 'NumChannels')
                        obj.NumChannels = varargin{i+1};
                    elseif strcmp(varargin{i}, 'Epsilon')
                        obj.Epsilon = double(varargin{i+1});
                    elseif strcmp(varargin{i}, 'Offset')
                        obj.Offset = double(varargin{i+1});
                    elseif strcmp(varargin{i}, 'Scale')
                        obj.Scale = double(varargin{i+1});
                    elseif strcmp(varargin{i}, 'NumInputs')
                        obj.NumInputs = double(varargin{i+1});
                    elseif strcmp(varargin{i}, 'InputNames')
                        obj.InputNames = varargin{i+1};
                    elseif strcmp(varargin{i}, 'NumOutputs')
                        obj.NumOutputs = double(varargin{i+1});
                    elseif strcmp(varargin{i}, 'OutputNames')
                        obj.OutputNames = varargin{i+1};
                    end
                    
                end
                if isempty(obj.Offset)
                    obj.Offset = zeros(1,obj.NumChannels);
                end
                if isempty(obj.Scale)
                    obj.Scale = ones(1,obj.NumChannels);
                end
                
            end
                
             
        end
        
        
    end
        
    % evaluation method
    methods
        
        function y = evaluate_old(obj, input)
            
        end 
        
        function y = evaluateSequence(obj, input)
            newInput = dlarray(input);
            if length(size(input))==2
                y = layernorm(newInput, obj.Offset, obj.Scale,"DataFormat",'CS','Epsilon',obj.Epsilon);
            elseif length(size(input))==3
                y = layernorm(newInput, obj.Offset, obj.Scale,"DataFormat",'SCS','Epsilon',obj.Epsilon);
            end
            y = extractdata(y);
        end
    end
    
    
        
    % exact reachability analysis using ImageStar or ImageZono
    methods
        
        
        function image = reach_star_single_input_old(obj, in_image)
            % @in_image: an input imagestar
            % @image: output set
            
            % author:Neelanjana Pal
            % date: 1/7/2020
            
            if ~isa(in_image, 'ImageStar')
                error('Input is not an ImageStar');
            end
            
            if isempty(obj.TrainedMean) || isempty(obj.TrainedVariance) || isempty(obj.Epsilon) || isempty(obj.Offset) || isempty(obj.Scale)
                error('Layer Normalization Layer does not have enough parameters');
            end
            
            var = obj.TrainedVariance;
            eps = obj.Epsilon;
            mean = obj.TrainedMean;
            scale = obj.Scale; 
            offset = obj.Offset;
            l(1,1, obj.NumChannels) = 0;
            for i=1:obj.NumChannels
                l(1,1,i) = 1/sqrt(var(1,1,i) + eps);
            end
            
            x = in_image.affineMap(l, -l.*mean);
            image = x.affineMap(scale, offset);
            
        end
        
        function image = reach_star_single_input(obj, in_image)
            % @in_image: an input imagestar
            % @image: output set
            
            % author: Mykhailo Ivashchenko
            % date: 9/17/2022
            
            if ~isa(in_image, 'ImageStar') && ~isa(in_image, 'Star') % CHANGED
                error('Input is not a Star or ImageStar');
            end
                       
            var = obj.TrainedVariance;
            eps = obj.Epsilon;
            mean = obj.TrainedMean;
            scale = obj.Scale; 
            offset = obj.Offset;
                        
            image = in_image;
            
            if isa(image, 'ImageStar')
                if(length(size(obj.TrainedMean)) == 2) && size(image.V, 1) == 1 && size(image.V, 2) == 1 && length(size(image.V)) == 4
                    obj.TrainedMean = reshape(obj.TrainedMean, [1 1 size(obj.TrainedMean, 1)]);
                    obj.TrainedVariance = reshape(obj.TrainedVariance, [1 1 size(obj.TrainedVariance, 1)]);
                    obj.Epsilon = reshape(obj.Epsilon, [1 1 size(obj.Epsilon, 1)]);
                    
                    if length(size(obj.Offset)) ~= 3
                        obj.Offset = reshape(obj.Offset, [1 1 size(obj.Offset, 1)]);
                    end
                    
                    if length(size(obj.Scale)) ~= 3
                        obj.Scale = reshape(obj.Scale, [1 1 size(obj.Scale, 1)]);
                    end
                    obj.NumChannels = 1;
                end
                
                x = in_image;
                
                if ~isempty(obj.TrainedMean) && ~isempty(obj.TrainedVariance) && ~isempty(obj.Epsilon) && ~isempty(obj.Offset) && ~isempty(obj.Scale)
                    l(1,1, obj.NumChannels) = 0;
                    for i=1:obj.NumChannels
                        l(1,1,i) = 1/sqrt(obj.TrainedVariance(1,1,i) + obj.Epsilon);
                    end
                    x = x.affineMap(l, -l.*obj.TrainedMean);
                end   
                
                image = x.affineMap(obj.Scale, obj.Offset);
            elseif isa(image, 'Star')
                l = scale ./ sqrt(var + eps);
                
                for i=1:size(image.V, 2)
                    image.V(:, i) = ((image.V(:, i) - mean) .* l + offset);
                end
            end
            
        end
        
        function image = reach_zono(obj, in_image)
            % @in_image: an input ImageZono
            % @image: output set
            
            % author:Neelanjana Pal
            %  date: 1/7/2020
            
            if ~isa(in_image, 'ImageZono')
                error('Input is not an ImageZono');
            end
            
            if isempty(obj.TrainedMean) || isempty(obj.TrainedVariance) || isempty(obj.Epsilon) || isempty(obj.Offset) || isempty(obj.Scale)
                error('Layer Normalization Layer does not have enough parameters');
            end
            
            var = obj.TrainedVariance;
            eps = obj.Epsilon;
            mean = obj.TrainedMean;
            scale = obj.Scale; 
            offset = obj.Offset;
            l(1,1, obj.NumChannels) = 0;
            for i=1:obj.NumChannels
                l(1,1,i) = 1/sqrt(var(1,1,i) + eps);
            end
            
            x = in_image.affineMap(l, -l.*mean);
            image = x.affineMap(scale, offset);
            
        end
        
        function images = reach_star_multipleInputs(obj, in_images, option)
            % @in_images: an array of ImageStars input set
            % @option: = 'parallel' or 'single' or empty
            
            % author:Neelanjana Pal
            % date: 1/7/2020
            
            n = length(in_images);
            if isa(in_images(n), 'ImageStar')
                images(n) = ImageStar; 
            else
                images(n) = Star; 
            end 
            
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

            % author:Neelanjana Pal
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

        function image = reach_star_single_input_Sequence(obj, in_image)
            % @in_image: an input imagestar
            % @image: output set
            
            % author: Neelanjana Pal
            % date: 1/17/2023
            
            if ~isa(in_image, 'ImageStar')
                error('The input is not an ImageStar object');
            end
            
            % compute output sets
            
            if isempty(in_image.im_lb) && isempty(in_image.im_ub)
                c = obj.evaluateSequence(in_image.V(:,:,:,1));
                layer = obj;
                layer.Offset = zeros(obj.NumChannels,1);
                parfor i = 2: size(in_image.V,4)
                    V(:,:,:,i) = layer.evaluateSequence(in_image.V(:,:,:,i));
                end
                V(:,:,:,1) = c;
                image = ImageStar(V, in_image.C, in_image.d, in_image.pred_lb, in_image.pred_ub);
            else
                im_lb = in_image.im_lb;
                im_ub = in_image.im_ub;
                lb = obj.evaluateSequence(im_lb);
                ub = obj.evaluateSequence(im_ub);
    
                image = ImageStar(lb,ub);
            end 
        end
        

        function images = reach_star_multipleInputs_Sequence(obj, in_images, option)
            % @in_images: an array of ImageStars input set
            % @option: = 'parallel' or 'single' or empty
            
            % author: Neelanjana Pal
            % date: 17/1/2023
            
            n = length(in_images);
            images(n) = ImageStar;  
            
            if strcmp(option, 'parallel')
                parfor i=1:n
                    images(i) = obj.reach_star_single_input_Sequence(in_images(i));
                end
            elseif strcmp(option, 'single') || isempty(option)
                for i=1:n
                    images(i) = obj.reach_star_single_input_Sequence(in_images(i));
                end
            else
                error('Unknown computation option');

            end
        end
       
        
        function images = reach(varargin)
            % @in_image: an input imagestar or imagezono
            % @image: output set
            % @option: = 'single' or 'parallel' 
            
            % author:Neelanjana Pal
            % date: 6/26/2019
             
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
                    option = [];
                otherwise
                    error('Invalid number of input arguments (should be 2, 3, 4, 5, or 6)');
            end
            
            
            if strcmp(method, 'approx-star') || strcmp(method, 'exact-star') || strcmp(method, 'abs-dom')|| contains(method, "relax-star")
                images = obj.reach_star_multipleInputs(in_images, option);
            elseif strcmp(method, 'approx-zono')
                images = obj.reach_zono_multipleInputs(in_images, option);
            else
                error("Unknown reachability method");
            end
          
        end
                 
        function images = reachSequence(varargin)
        %obj = varargin{1};
        %seqs = obj.reach(varargin{2:nargin});
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
                    option = [];
                otherwise
                    error('Invalid number of input arguments (should be 2, 3, 4, 5, or 6)');
            end

            if strcmp(method, 'approx-star') || strcmp(method, 'exact-star') || strcmp(method, 'abs-dom')|| contains(method, "relax-star")
                images = obj.reach_star_multipleInputs_Sequence(in_images, option);
            elseif strcmp(method, 'approx-zono')
                images = obj.reach_zono_multipleInputs(in_images, option); %need to change
            else
                error("Unknown reachability method");
            end
            

        end
    end
    
    
    methods(Static)
         % parse a trained batch normalization layer from matlab
        function L = parse(layer)
            % @layer: batch normalization layer
            % @L: constructed layer
                        
            % author:Neelanjana Pal
            % date: 8/30/2023
            
            
            if ~isa(layer, 'nnet.cnn.layer.LayerNormalizationLayer')
                error('Input is not a Matlab nnet.cnn.layer.LayerNormalizationLayer class');
            end
                       
            L = LayerNormalizationLayer('Name', layer.Name, 'NumChannels', layer.NumChannels,'Epsilon', layer.Epsilon, 'Offset', layer.Offset, 'Scale', layer.Scale, 'NumInputs', layer.NumInputs, 'InputNames', layer.InputNames, 'NumOutputs', layer.NumOutputs, 'OutputNames', layer.OutputNames);
            
        end
        
    end
    
    
    
end

