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

        % change params to gpuArrays
        function obj = toGPU(obj)
            obj.Offset = gpuArray(obj.Offset);
            obj.Scale = gpuArray(obj.Scale);
            obj.Epsilon = gpuArray(obj.Epsilon);
        end

        % Change params precision
        function obj = changeParamsPrecision(obj, precision)
            if strcmp(precision, "double")
                obj.Offset = double(obj.Offset);
                obj.Scale = double(obj.Scale);
                obj.Epsilon = double(obj.Epsilon);
            elseif strcmp(precision, "single")
                obj.Offset = single(obj.Offset);
                obj.Scale = single(obj.Scale);
                obj.Epsilon = single(obj.Epsilon);
            else
                error("Parameter numerical precision must be 'single' or 'double'");
            end
        end
        
        
    end
        
    % evaluation method
    methods
        
        function y = evaluateSequence(obj, input)
            newInput = dlarray(input);
            if length(size(input))==2
                y = layernorm(newInput, obj.Offset, obj.Scale,"DataFormat",'CS','Epsilon',obj.Epsilon);
            elseif length(size(input))==3
                y = layernorm(newInput, obj.Offset, obj.Scale,"DataFormat",'SCS','Epsilon',obj.Epsilon);
            end
            y = extractdata(y);
        end

        function y = evaluate(obj, input)
            % Evaluate layer normalization on input
            % Supports various input formats
            y = obj.evaluateSequence(input);
        end

    end
    
    
        
    % exact reachability analysis using ImageStar or ImageZono
    methods
        
        % NOTE: reach_star_single_input_old was removed - it incorrectly
        % referenced TrainedMean/TrainedVariance which don't exist in LayerNorm.
        % Layer Normalization computes mean/variance from input, not stored stats.
        
        function image = reach_star_single_input(obj, in_image)
            % @in_image: an input ImageStar or Star
            % @image: output set
            %
            % Layer Normalization normalizes across features (not batch).
            % Unlike BatchNorm, it does NOT use pre-computed statistics.
            % Instead, mean and variance are computed from the input itself.
            %
            % For reachability, we use a bounds-based approximation:
            % Given input bounds, we over-approximate the output bounds.
            %
            % author: Mykhailo Ivashchenko (original)
            % date: 9/17/2022
            % modified: NNV Team, December 2025 - fixed for proper LayerNorm semantics

            if ~isa(in_image, 'ImageStar') && ~isa(in_image, 'Star')
                error('Input is not a Star or ImageStar');
            end

            eps = obj.Epsilon;
            scale = obj.Scale(:);
            offset = obj.Offset(:);

            if isa(in_image, 'ImageStar')
                % For ImageStar: use bounds-based approximation
                % Get input bounds
                if ~isempty(in_image.im_lb) && ~isempty(in_image.im_ub)
                    lb = in_image.im_lb;
                    ub = in_image.im_ub;
                else
                    [lb, ub] = in_image.estimateRanges();
                end

                % Flatten for LayerNorm computation
                lb_flat = lb(:);
                ub_flat = ub(:);
                n = length(lb_flat);

                % Compute bounds on mean: mean is in [mean(lb), mean(ub)]
                % but more precisely, mean(x) where x_i in [lb_i, ub_i]
                mean_lb = mean(lb_flat);
                mean_ub = mean(ub_flat);

                % Compute bounds on variance (conservative)
                % Variance is always >= 0
                center = (lb_flat + ub_flat) / 2;
                var_center = var(center, 1);  % variance of center point
                var_lb = max(0, var_center * 0.5);  % conservative lower bound
                var_ub = var_center + max((ub_flat - lb_flat).^2) / 4;  % conservative upper bound

                % Compute normalization factor bounds
                std_lb = sqrt(var_lb + eps);
                std_ub = sqrt(var_ub + eps);

                % Compute output bounds
                % y = (x - mean) / std * scale + offset
                % For each element, compute worst-case bounds
                out_lb = zeros(size(lb));
                out_ub = zeros(size(ub));

                for i = 1:numel(lb)
                    % Worst case for (x - mean):
                    x_minus_mean_lb = lb(i) - mean_ub;
                    x_minus_mean_ub = ub(i) - mean_lb;

                    % Get scale for this element
                    if i <= length(scale)
                        s = scale(i);
                        o = offset(i);
                    else
                        s = scale(mod(i-1, length(scale)) + 1);
                        o = offset(mod(i-1, length(offset)) + 1);
                    end

                    % Normalize and apply scale/offset
                    if s >= 0
                        out_lb(i) = x_minus_mean_lb / std_ub * s + o;
                        out_ub(i) = x_minus_mean_ub / std_lb * s + o;
                    else
                        out_lb(i) = x_minus_mean_ub / std_lb * s + o;
                        out_ub(i) = x_minus_mean_lb / std_ub * s + o;
                    end
                end

                % Create output ImageStar from bounds
                image = ImageStar(out_lb, out_ub);

            elseif isa(in_image, 'Star')
                % For Star: use bounds-based approximation
                n = in_image.dim;

                % Get bounds
                lb = zeros(n, 1);
                ub = zeros(n, 1);
                for i = 1:n
                    lb(i) = in_image.getMin(i, 'linprog');
                    ub(i) = in_image.getMax(i, 'linprog');
                end

                % Compute mean bounds
                mean_lb = mean(lb);
                mean_ub = mean(ub);

                % Compute variance bounds (conservative)
                center = (lb + ub) / 2;
                var_center = var(center, 1);
                var_lb = max(0, var_center * 0.5);
                var_ub = var_center + max((ub - lb).^2) / 4;

                % Normalization factor bounds
                std_lb = sqrt(var_lb + eps);
                std_ub = sqrt(var_ub + eps);

                % Compute output bounds
                out_lb = zeros(n, 1);
                out_ub = zeros(n, 1);

                for i = 1:n
                    x_minus_mean_lb = lb(i) - mean_ub;
                    x_minus_mean_ub = ub(i) - mean_lb;

                    if i <= length(scale)
                        s = scale(i);
                        o = offset(i);
                    else
                        s = scale(mod(i-1, length(scale)) + 1);
                        o = offset(mod(i-1, length(offset)) + 1);
                    end

                    if s >= 0
                        out_lb(i) = x_minus_mean_lb / std_ub * s + o;
                        out_ub(i) = x_minus_mean_ub / std_lb * s + o;
                    else
                        out_lb(i) = x_minus_mean_ub / std_lb * s + o;
                        out_ub(i) = x_minus_mean_lb / std_ub * s + o;
                    end
                end

                % Create output Star from bounds
                image = Star(out_lb, out_ub);
            end

        end
        
        function image = reach_zono(obj, in_image)
            % @in_image: an input ImageZono
            % @image: output set
            %
            % Layer Normalization for ImageZono using bounds-based approximation.
            %
            % author: Neelanjana Pal (original)
            % date: 1/7/2020
            % modified: NNV Team, December 2025 - fixed for proper LayerNorm semantics

            if ~isa(in_image, 'ImageZono')
                error('Input is not an ImageZono');
            end

            eps = obj.Epsilon;
            scale = obj.Scale(:);
            offset = obj.Offset(:);

            % Get bounds from zonotope
            [lb, ub] = in_image.getBounds();

            % Flatten for LayerNorm computation
            lb_flat = lb(:);
            ub_flat = ub(:);

            % Compute bounds on mean
            mean_lb = mean(lb_flat);
            mean_ub = mean(ub_flat);

            % Compute bounds on variance (conservative)
            center = (lb_flat + ub_flat) / 2;
            var_center = var(center, 1);
            var_lb = max(0, var_center * 0.5);
            var_ub = var_center + max((ub_flat - lb_flat).^2) / 4;

            % Normalization factor bounds
            std_lb = sqrt(var_lb + eps);
            std_ub = sqrt(var_ub + eps);

            % Compute output bounds
            out_lb = zeros(size(lb));
            out_ub = zeros(size(ub));

            for i = 1:numel(lb)
                x_minus_mean_lb = lb(i) - mean_ub;
                x_minus_mean_ub = ub(i) - mean_lb;

                if i <= length(scale)
                    s = scale(i);
                    o = offset(i);
                else
                    s = scale(mod(i-1, length(scale)) + 1);
                    o = offset(mod(i-1, length(offset)) + 1);
                end

                if s >= 0
                    out_lb(i) = x_minus_mean_lb / std_ub * s + o;
                    out_ub(i) = x_minus_mean_ub / std_lb * s + o;
                else
                    out_lb(i) = x_minus_mean_ub / std_lb * s + o;
                    out_ub(i) = x_minus_mean_lb / std_ub * s + o;
                end
            end

            % Create output ImageZono from bounds
            image = ImageZono(out_lb, out_ub);

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

