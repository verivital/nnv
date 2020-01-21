classdef ImageZono < handle
    % Class for representing set of images using Zonotope
    % An image can be attacked by bounded noise. An attacked image can
    % be represented using an ImageZono Set
    % Dung Tran: 1/2/2020
   
    %=================================================================%
    
    
    properties
        numChannels = 0; % number of channels, e.g., color images have 3 channel
        height = 0; % height of image
        width = 0; % width of image
        
        % A box representation of an ImageZono
        % A convenient way for user to specify the attack
        
        lb_image = []; % lower bound of attack (high-dimensional array)
        ub_image = []; % upper bound of attack (high-dimensional array)
        
        %                      2-Dimensional ImageZono
        % ====================================================================%
        %                   Definition of 2-Dimensonal ImageZono
        % 
        % A ImageZono Z= <c, V> is defined by: 
        % S = {x| x = c + a[1]*V[1] + a[2]*V[2] + ... + a[n]*V[n]
        %           = V * b, V = {c V[1] V[2] ... V[n]}, 
        %                    b = [1 a[1] a[2] ... a[n]]^T                                   
        %                    where -1 <= a[i] <= 1}
        % where, V[0], V[i] are 2D matrices with the same dimension, i.e., 
        % V[i] \in R^{m x n}
        % V[0] : is called the center matrix and V[i] is called the basic matrix 
        % [a[1]...a[n] are called predicate variables
        % The notion of 2D ImageZono is more general than the original Zonotope where
        % the V[0] and V[i] are vectors. 
        % 
        % Dimension of 2D ImageZono is the dimension of the center matrix V[0]
        % 
        % ====================================================================%
        % The 2D representation of ImageZono is convenient for reachability analysis
        V = []; % a cell (size = numPred)
        numPreds = 0; % number of predicate variables
        
    end
    
    methods
        % constructor using 2D representation of an ImageZono
        function obj = ImageZono(varargin)
            % @V:  an array of basis images          
            % @LB: lower bound of attack (high-dimensional array)
            % @UB: upper bound of attack (high-dimensional array)
                        
            % author: Dung Tran
            % date: 1/3/2019
            
            switch nargin
                case 2
                    obj.lb_image = varargin{1};
                    obj.ub_image = varargin{2};
                    
                    % check consistency
                    size1 = size(obj.lb_image);
                    size2 = size(obj.ub_image);

                    cp12 = (size1 == size2);

                    if sum(cp12) ~= 3 && sum(cp12) ~= 2
                        error('Different sizes between lower bound image and upper bound image');
                    end

                    obj.height = size1(1);
                    obj.width = size1(2);
                    if length(size1) == 2
                        obj.numChannels = 1;
                    else
                        obj.numChannels = size1(3);
                    end
                    
                    % get basis images array
                    
                    lb = reshape(obj.lb_image, [obj.height*obj.width*obj.numChannels 1]);
                    ub = reshape(obj.ub_image, [obj.height*obj.width*obj.numChannels 1]);
                    
                    S = Star(lb, ub);
                    obj.numPreds = S.nVar;
                    obj.V = reshape(S.V, [obj.height obj.width obj.numChannels obj.numPreds + 1]);
                    
                    
                case 1
                    
                    obj.V = varargin{1};
                    [obj.height, obj.width, obj.numChannels] = size(obj.V(:,:,:,1));
                    obj.numPreds = size(obj.V, 4) - 1;
                    
                    center = obj.V(:,:,:,1);
                    generators = obj.V(:,:,:, 2:obj.numPreds + 1);
                    center = reshape(center, [obj.height*obj.width*obj.numChannels 1]);
                    generators = reshape(generators, [obj.height*obj.width*obj.numChannels obj.numPreds]);

                    Z = Zono(center, generators);
                    [lb, ub] = Z.getBounds;
                    
                    obj.lb_image = reshape(lb, [obj.height obj.width obj.numChannels]);
                    obj.ub_image = reshape(ub, [obj.height obj.width obj.numChannels]);
                    
                case 0
                    
                    
                    
                otherwise
                    error('Invalid number of inputs, shoule be 0, 1 or 2');
            end
            
            
            
            
   
        end
                
        % randomly generate a set of images from an imageZono set
        function images = sample(obj, N)
            % @N: number of images 
            
            % author: Dung Tran
            % date: 1/3/2020
            
            if isempty(obj.V)
                error('The imagezono is an empty set');
            end
            
            % fill code here
            images = {};
              
        end
        
        
        % evaluate an ImageZono with specific values of predicates
        function image = evaluate(obj, pred_val)
            % @pred_val: valued vector of predicate variables
            
            % author: Dung Tran
            % date: 1/3/2020
            
            if isempty(obj.V)
                error('The ImageZono is an empty set');
            end
            
            if size(pred_val, 2) ~= 1
                error('Invalid predicate vector');
            end
            
            if size(pred_val, 1) ~= obj.numPreds
                error('Inconsistency between the size of the predicate vector and the number of predicates in the ImageZono');
            end
            
            % check if all values of predicate variables are in [-1, 1]           
            for i=1:obj.numPreds
                if ~(pred_val(i)<= 1 && pred_val(i) >= -1)
                    error('Predicate values should be in the range of [-1, 1] for ImageZono');
                end
            end
            
            
            image(:, :, obj.numChannels) = zeros(obj.height, obj.width);
            for i=1:obj.numChannels
                image(:, :, i) = obj.V(:,:,i, 1);
                for j=2:obj.numPreds + 1
                    image(:, :, i) = image(:, :, i) + pred_val(j-1) * obj.V(:,:,i, j);
                end
            end
                      
        end
        
        
                
        % affineMap of an ImageZono is another imagezono
        % y = scale * x + offset;
        function image = affineMap(obj, scale, offset)
            % @scale: scale coefficient [1 x 1 x NumChannels] array
            % @offset: offset coefficient [1 x 1 x NumChannels] array
            % @image: a new ImageZono
            
            % author: Dung Tran
            % date: 1/1/2020
            
            
            if ~isempty(scale) && ~isscalar(scale) && size(scale, 3) ~= obj.numChannels
                error('Inconsistent number of channels between scale array and the ImageStar');
            end
                        
            
            if ~isempty(scale) 
                new_V = scale.*obj.V;
            else
                new_V = obj.V;
            end
            
            if ~isempty(offset)
                new_V(:,:,:,1) = new_V(:,:,:,1) + offset;
            end
                       
            image = ImageZono(new_V);
                     
        end
        
        % transform to Zono
        function Z = toZono(obj)
            
            % author: Dung Tran
            % date: 1/2/2020
            
            center = obj.V(:,:,:,1);
            generators = obj.V(:,:,:, 2:obj.numPreds + 1);
            
            center = reshape(center, [obj.height*obj.width*obj.numChannels 1]);
            generators = reshape(generators, [obj.height*obj.width*obj.numChannels obj.numPreds]);
            
            Z = Zono(center, generators);
            
        end
        
        % transform to ImageStar
        function S = toImageStar(obj)
            % author: Dung Tran
            % date: 1/2/2020
            
            pred_lb = -ones(obj.numPreds, 1);
            pred_ub = ones(obj.numPreds, 1); 
            
            C1 = eye(obj.numPreds);
            d1 = pred_ub;
            C2 = -C1;
            d2 = pred_ub; 
            
            C = [C1; C2];
            d = [d1; d2]; 
            
            S = ImageStar(obj.V, C, d, pred_lb, pred_ub, obj.lb_image, obj.ub_image);
            
        end
        
        
        % contain, check if an ImageZono contain an image 
        function bool = contains(obj, image)
            % @image: input image
            % @bool: = 1 if the ImageStar contain the image
            %        = 0 if the ImageStar does not contain the image
            
            % author: Dung Tran
            % date: 1/8/2020
            
            n = size(image);
            if length(n) == 2 % one channel image
                if n(1) ~= obj.height || n(2) ~= obj.width || obj.numChannels ~= 1
                    error('Inconsistent dimenion between input image and the ImageStar');
                end
                y = reshape(image, [n(1)*n(2) 1]);
            elseif length(n) == 3
                if n(1) ~= obj.height || n(2) ~= obj.width || n(3) ~= obj.numChannels
                    error('Inconsistent dimenion between input image and the ImageStar');
                end
                y = reshape(image, [n(1)*n(2)*n(3) 1]);
            else
                error('Invalid input image');
            end
            
            Z = obj.toZono;
            bool = Z.contains(y);
            
        end
        
        % get Ranges
        function [lb, ub] = getRanges(obj)
            % author: Dung Tran
            % date: 1/9/2020
            
            lb = obj.lb_image;
            ub = obj.ub_image;
        end
        
        function b = is_p1_larger_p2(obj, p1, p2)
            % @p1: the first point = [h1, w1, c1]
            % @p2: the second point = [h2, w2, c2]
            % h: height, w: width, c: channel index
            
            % @b = 1 -> p1 > p2 is feasible
            %    = 0 -> p1 > p2 is not feasible
            
            % author: Dung Tran
            % date: 7/23/2019
            
            % a*(V2 - V1) <= c1 - c2
            
            S = obj.toImageStar;
            b = S.is_p1_larger_p2(p1, p2);           
            
        end
        
        
       
         
    end
    
    
    
    methods(Static)
        
        
    end
    
    
    
    
    
    
   
        
       
end

