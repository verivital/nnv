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

                % SOUND output bounds. The previous code used var_center*0.5 /
                % max(.)^2/4 -- UNSOUND: it took the variance of the CENTER
                % point (0 even when the true variance is large), giving
                % std_lb=sqrt(eps) and a lower bound that EXCLUDED reachable
                % outputs (2532/5000 MC violations; counterexample x=[0;5] over
                % [0,2]x[-3,5]). See sound_bounds + LayerNorm soundness test.
                [olb, oub] = obj.sound_bounds(lb_flat, ub_flat, scale, offset, eps);
                image = ImageStar(reshape(olb, size(lb)), reshape(oub, size(ub)));

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

                % SOUND output bounds (see sound_bounds; the old var_center
                % formulas were unsound -- same bug as the ImageStar path).
                [olb, oub] = obj.sound_bounds(lb, ub, scale, offset, eps);
                image = Star(olb, oub);
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

            % SOUND output bounds (same sound bound as the Star path; the old
            % var_center formula was unsound here too).
            [olb, oub] = obj.sound_bounds(lb_flat, ub_flat, scale, offset, eps);
            image = ImageZono(reshape(olb, size(lb)), reshape(oub, size(ub)));

        end

        function [out_lb, out_ub] = sound_bounds(~, lb, ub, scale, offset, eps)
            % SOUND interval over-approximation of LayerNorm output bounds.
            %   y_i = scale_i * (x_i - mean)/sqrt(var + eps) + offset_i,
            %   mean = (1/n) sum_j x_j,  var = (1/n) sum_j (x_j - mean)^2.
            % Combines two sound over-approximations of the normalized z_i and
            % intersects them:
            %  (1) Input-dependent interval: mean in [mean(lb),mean(ub)]; a sound
            %      UPPER bound on var via the per-coordinate worst-case squared
            %      deviation; var lower bound 0 (=> std_lo = sqrt(eps)); then
            %      z_i = (x_i-mean)/std by interval division with num and std
            %      taken independently (a sound SUPERSET of the coupled range).
            %  (2) Universal cap: sum_j z_j = 0 and sum_j z_j^2 <= n imply
            %      |z_i| <= sqrt(n-1) for ANY input and any eps >= 0.
            % (1) can blow up (std_lo ~ sqrt(eps)); intersecting with (2) keeps
            % the result both SOUND and finite. n=1 -> zcap=0 -> output=offset.
            lb = double(lb(:)); ub = double(ub(:));
            n = numel(lb);
            scale = scale(:); offset = offset(:);
            if numel(scale)  ~= n, scale  = scale(mod((0:n-1), numel(scale))  + 1); end
            if numel(offset) ~= n, offset = offset(mod((0:n-1), numel(offset)) + 1); end

            mean_lb = mean(lb); mean_ub = mean(ub);
            d_lo = lb - mean_ub; d_hi = ub - mean_lb;       % x_j - mean range
            var_ub = mean(max(d_lo.^2, d_hi.^2));            % sound upper bound on var
            std_lo = sqrt(eps);                              % var_lb = 0 (sound)
            std_hi = sqrt(var_ub + eps);
            zcap = sqrt(max(n - 1, 0));                      % universal |z_i| bound

            num_lo = lb - mean_ub; num_hi = ub - mean_lb;    % x_i - mean
            z_lo = zeros(n,1); z_hi = zeros(n,1);
            for i = 1:n
                a = num_lo(i); b = num_hi(i);
                if a >= 0
                    z_lo(i) = a / std_hi;  z_hi(i) = b / std_lo;
                elseif b <= 0
                    z_lo(i) = a / std_lo;  z_hi(i) = b / std_hi;
                else
                    z_lo(i) = a / std_lo;  z_hi(i) = b / std_lo;
                end
            end
            z_lo = max(z_lo, -zcap);                         % intersect with cap
            z_hi = min(z_hi,  zcap);

            out_lb = zeros(n,1); out_ub = zeros(n,1);
            pos = scale >= 0;
            out_lb(pos)  = scale(pos).*z_lo(pos)   + offset(pos);
            out_ub(pos)  = scale(pos).*z_hi(pos)   + offset(pos);
            out_lb(~pos) = scale(~pos).*z_hi(~pos) + offset(~pos);
            out_ub(~pos) = scale(~pos).*z_lo(~pos) + offset(~pos);
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
            
            % SOUND output set. The previous code was UNSOUND: it applied the
            % NONLINEAR LayerNorm to im_lb/im_ub (or to each basis vector V(:,:,:,i))
            % independently, as if LayerNorm were affine/monotone -- but LayerNorm
            % couples elements and is non-monotone, so f(lb)/f(ub) and the
            % per-generator image do NOT bound the true reachable set.
            % LayerNorm normalizes over the channel dim, so the normalized value
            % z_i = (x_i - mean)/sqrt(var + eps) obeys the UNIVERSAL bound
            % |z_i| <= sqrt(C-1) (C = NumChannels = group size) for ANY input and
            % any position -> y_i in offset_i +/- |scale_i|*sqrt(C-1). Sound (input-
            % independent; loose). A tighter input-dependent per-group bound is
            % future work (see sound_bounds + the LayerNorm grouping note).
            C = obj.NumChannels;
            zcap = sqrt(max(C - 1, 0));
            scale = obj.Scale(:); offset = obj.Offset(:);
            if ~isempty(in_image.im_lb)
                sz = size(in_image.im_lb);
            else
                sz = size(in_image.V(:,:,:,1));
            end
            ne = prod(sz);
            sc = scale(mod((0:ne-1).', numel(scale))  + 1);
            of = offset(mod((0:ne-1).', numel(offset)) + 1);
            lo = of - abs(sc) * zcap;
            hi = of + abs(sc) * zcap;
            image = ImageStar(reshape(lo, sz), reshape(hi, sz));
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

