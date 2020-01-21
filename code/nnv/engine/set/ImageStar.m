classdef ImageStar < handle
    % Class for representing set of images using Star set
    % An image can be attacked by bounded noise. An attacked image can
    % be represented using an ImageStar Set
    % Dung Tran: 6/12/2019
    
    %=================================================================%
    %   a 3-channels color image is represented by 3-dimensional array 
    %   Each dimension contains a h x w matrix, h and w is the height
    %   width of the image. h * w = number of pixels in the image.
    %   *** A gray image has only one channel.
    %
    %   Problem: How to represent a disturbed(attacked) image?
    %   
    %   Use a center image (a matrix) + a disturbance matrix (positions
    %   of attacks and bounds of corresponding noises)
    %
    %   For example: Consider a 4 x 4 (16 pixels) gray image 
    %   The image is represented by 4 x 4 matrix:
    %               IM = [1 1 0 1; 0 1 0 0; 1 0 1 0; 0 1 1 1]
    %   This image is attacked at pixel (1,1) (1,2) and (2,4) by bounded
    %   noises:     |n1| <= 0.1, |n2| <= 0.2, |n3| <= 0.05
    %
    %
    %   Lower and upper noises bounds matrices are: 
    %         LB = [-0.1 -0.2 0 0; 0 0 0 -0.05; 0 0 0 0; 0 0 0 0]
    %         UB = [0.1 0.2 0 0; 0 0 0 0.05; 0 0 0 0; 0 0 0 0]
    %   The lower and upper bounds matrices also describe the position of 
    %   attack.
    %
    %   Under attack we have: -0.1 + 1 <= IM(1,1) <= 1 + 0.1
    %                         -0.2 + 1 <= IM(1,2) <= 1 + 0.2
    %                            -0.05 <= IM(2,4) <= 0.05
    %
    %   To represent the attacked image we use IM, LB, UB matrices
    %   For multi-channel image we use multi-dimensional array IM, LB, UB
    %   to represent the attacked image. 
    %   For example, for an attacked color image with 3 channels we have
    %   IM(:, :, 1) = IM1, IM(:,:,2) = IM2, IM(:,:,3) = IM3
    %   LB(:, :, 1) = LB1, LB(:,:,2) = LB2, LB(:,:,3) = LB3
    %   UB(:, :, 1) = UB1, UB(:,:,2) = UB2, UB(:,:,3) = UB3
    %   
    %   The image object is: image = ImageStar(IM, LB, UB)
    %=================================================================%
    
    
    properties
        numChannel = 0; % number of channels, e.g., color images have 3 channel
        height = 0; % height of image
        width = 0; % width of image
        
        % A box representation of an ImageStar
        % A convenient way for user to specify the attack
        
        IM = []; % center image (high-dimensional array)
        LB = []; % lower bound of attack (high-dimensional array)
        UB = []; % upper bound of attack (high-dimensional array)
        
        % 2D representation of an ImageStar
        % ====================================================================%
        %                   Definition of Star2D
        % 
        % A 2D star set S is defined by: 
        % S = {x| x = V[0] + a[1]*V[1] + a[2]*V[2] + ... + a[n]*V[n]
        %           = V * b, V = {c V[1] V[2] ... V[n]}, 
        %                    b = [1 a[1] a[2] ... a[n]]^T                                   
        %                    where C*a <= d, constraints on a[i]}
        % where, V[0], V[i] are 2D matrices with the same dimension, i.e., 
        % V[i] \in R^{m x n}
        % V[0] : is called the center matrix and V[i] is called the basic matrix 
        % [a[1]...a[n] are called predicate variables
        % C: is the predicate constraint matrix
        % d: is the predicate constraint vector
        %
        % The notion of Star2D is more general than the original Star set where
        % the V[0] and V[i] are vectors. 
        % 
        % Dimension of Star2D is the dimension of the center matrix V[0]
        % 
        % ====================================================================%
        % The 2D representation of ImageStar is convenient for reachability analysis
        V = []; % a cell (size = numPred)
        C = []; % a constraints matrix of the predicate
        d = []; % a constraints vector of the predicate
        numPred = 0; % number of predicate variables
        pred_lb = []; % lower bound vector of the predicate
        pred_ub = []; % upper bound vector of the predicate
        im_lb = []; % lower bound image of the ImageStar
        im_ub = []; % upper bound image of the ImageStar   

    end
    
    methods
        % constructor using 2D representation/1D representation of an ImageStar
        function obj = ImageStar(varargin)
            % @nargin = 3: IM = varargin{1}, LB = varagin{2}, UB =
            % varargin{3}
            %         = 2: Stars = varargin{1}, imageSize = varagin{2}
            %         = otherwise: IM = [], LB = [], UB = [], Stars = []
            % @IM: center image (high-dimensional array)
            % @LB: lower bound of attack (high-dimensional array)
            % @UB: upper bound of attack (high-dimensional array)
            % @Stars: 1D representation of an ImageStar (flattened ImageStar)
            
            % author: Dung Tran
            % date: 12/17/2018
            
            switch nargin
                
                case 3 % input center image and lower and upper bound matrices (box-representation)
                    
                    IM1 = varargin{1};
                    LB1 = varargin{2};
                    UB1 = varargin{3};
                    n = size(IM1); % n(1) and n(2) are height and width of image
                                  % n(3) is number of channels
                    l = size(LB1);
                    u = size(UB1);
                    
                    if n(1) ~= l(1) || n(1) ~= u(1) || n(2) ~= l(2) || n(2) ~= u(2) 
                        error('Inconsistency between center image and attack bound matrices');
                    end

                    if length(n) ~= length(l) || length(n) ~= length(u)
                        error('Inconsistency between center image and attack bound matrices');
                    end

                    if length(n) == 2 && length(l) == 2 && length(u) == 2

                        obj.numChannel = 1;
                        obj.IM = IM1;
                        obj.LB = LB1;
                        obj.UB = UB1;
                        obj.height = n(1);
                        obj.width = n(2);

                    elseif length(n) == 3 && length(l) == 3 && length(u) == 3

                        if n(3) == l(3) && n(3) == u(3)
                            obj.numChannel = n(3);
                            obj.IM = IM1; 
                            obj.LB = LB1;
                            obj.UB = UB1;
                            obj.height = n(1);
                            obj.width = n(2);
                        else
                            error('Inconsistent number of channels between the center image and the bound matrices');
                        end

                    else
                        error('Inconsistent number of channels between the center image and the bound matrices');
                    end
                    
                    obj.im_lb = IM1 + LB1; % lower bound image
                    obj.im_ub = IM1 + UB1; % upper bound image

                    % converting box ImageStar to an array of 2D Stars
                    
                    n = size(obj.im_lb);
                    I = Star(reshape(obj.im_lb, [n(1)*n(2)*n(3), 1]),reshape(obj.im_ub, [n(1)*n(2)*n(3), 1]));
                    
                    obj.V = reshape(I.V,[n(1), n(2), n(3), I.nVar + 1]);
                    obj.C = I.C;
                    obj.d = I.d;
                    obj.pred_lb = I.predicate_lb;
                    obj.pred_ub = I.predicate_ub;
                    obj.numPred = I.nVar;
 
                                         
                case 5 % 
                    
                    V1 = varargin{1};  % basis matrices
                    C1 = varargin{2};  % predicate constraint matrix 
                    d1 = varargin{3};  % predicate constraint vector
                    lb1 = varargin{4}; % predicate lower bound
                    ub1 = varargin{5}; % predicate upper bound
                    
                                        
                    if size(C1, 1) ~= size(d1, 1)
                        error('Inconsistent dimension between constraint matrix and constraint vector');
                    end
                    
                    if size(d1, 2) ~= 1
                        error('Invalid constraint vector, vector should have one column');
                    end
                                        
                    obj.numPred = size(C1, 2);
                    obj.C = C1;
                    obj.d = d1; 
                    
                    if size(C1, 2) ~= size(lb1, 1) || size(C1, 2) ~= size(ub1, 1)
                        error('Number of predicates is different from the size of the lower bound or upper bound predicate vector');
                    end
                    
                    if size(lb1, 2) ~= 1 || size(ub1, 2) ~= 1
                        error('Invalid lower/upper bound predicate vector, vector should have one column');
                    end
                    
                    obj.pred_lb = lb1;
                    obj.pred_ub = ub1;
                    
                    n = size(V1);
                    
                    if length(n) == 3
                        obj.V = V1;
                        obj.height = n(1);
                        obj.width = n(2);
                        obj.numChannel = n(3);
                        
                    elseif length(n) == 4
                        
                        if n(4) ~= obj.numPred + 1
                            error('Inconsistency between the basis matrix and the number of predicate variables');
                        else
                            obj.numChannel = n(3);
                            obj.V = V1;
                            obj.height = n(1);
                            obj.width = n(2);
                        end
                        
                    elseif length(n) == 2
                            obj.numChannel = 1;
                            obj.V = V1;
                            obj.height = n(1);
                            obj.width = n(2);
                    else
                        error('Invalid basis matrix');
                    end
                                            
                     
                case 7 % 
                    
                    V1 = varargin{1};  % basis matrices
                    C1 = varargin{2};  % predicate constraint matrix 
                    d1 = varargin{3};  % predicate constraint vector
                    lb1 = varargin{4}; % predicate lower bound
                    ub1 = varargin{5}; % predicate upper bound
                    im_lb1 = varargin{6}; % lower bound image
                    im_ub1 = varargin{7}; % upper bound image
                                                          
                                        
                    if size(C1, 1) ~= size(d1, 1)
                        error('Inconsistent dimension between constraint matrix and constraint vector');
                    end
                    
                    if size(d1, 2) ~= 1
                        error('Invalid constraint vector, vector should have one column');
                    end
                                        
                    obj.numPred = size(C1, 2);
                    obj.C = C1;
                    obj.d = d1; 
                    
                    if size(C1, 2) ~= size(lb1, 1) || size(C1, 2) ~= size(ub1, 1)
                        error('Number of predicates is different from the size of the lower bound or upper bound predicate vector');
                    end
                    
                    if size(lb1, 2) ~= 1 || size(ub1, 2) ~= 1
                        error('Invalid lower/upper bound predicate vector, vector should have one column');
                    end
                    
                    obj.pred_lb = lb1;
                    obj.pred_ub = ub1;
                 
                    n = size(V1);
                    
                    if length(n) == 3
                        
                        obj.numPred = 0;
                        obj.V = V1;
                        obj.height = n(1);
                        obj.width = n(2);
                        obj.numChannel = n(3);
                        
                    elseif length(n) == 4
                        
                        if n(4) ~= obj.numPred + 1
                            error('Inconsistency between the basis matrix and the number of predicate variables');
                        else
                            obj.numChannel = n(3);
                            obj.V = V1;
                            obj.height = n(1);
                            obj.width = n(2);
                        end
                        
                    elseif length(n) == 2
                            obj.numChannel = 1;
                            obj.numPred = 0;
                            obj.V = V1;
                            obj.height = n(1);
                            obj.width = n(2);
                    else
                        error('Invalid basis matrix');
                    end   
                    
                    if ~isempty(im_lb1) && (size(im_lb1,1) ~= obj.height || size(im_lb1, 2) ~= obj.width)
                        error('Inconsistent dimension between lower bound image and the constructed imagestar');
                    else
                        obj.im_lb = im_lb1;
                    end
                    
                    if ~isempty(im_ub1) && (size(im_ub1,1) ~= obj.height || size(im_ub1, 2) ~= obj.width)
                        error('Inconsistent dimension between upper bound image and the constructed imagestar');
                    else
                        obj.im_ub = im_ub1;
                    end
                    
                                                                    
                case 0 % create an empty ImageStar

                    obj.numChannel = 0; 
                    obj.height = 0;
                    obj.width = 0;
                    obj.IM = [];
                    obj.LB = [];
                    obj.UB = [];
                    obj.V = [];
                    obj.C = [];
                    obj.d = [];
                    obj.numPred = 0;
                    obj.pred_lb = [];
                    obj.pred_ub = [];
                    obj.im_lb = [];
                    obj.im_ub = [];
                                      
                otherwise
                    
                    error('Invalid number of input arguments, (should be from 0, 3, 5 , or 7)');
                    
            end
                        
        end
                
        % randomly generate a set of images from an imagestar set
        function images = sample(obj, N)
            % @N: number of images 
            
            % author: Dung Tran
            % date: 6/14/2019
            
            if isempty(obj.V)
                error('The imagestar is an empty set');
            end
            
            if isempty(obj.C) || isempty(obj.d)
                images = obj.IM;
            else
                
                V1 = [zeros(obj.numPred, 1) eye(obj.numPred)];
                S = Star(V1, obj.C, obj.d);                
                pred_samples = S.sample(N); 
                
                M = length(pred_samples);
                images = cell(1, M);
                for i=1:M
                    images{i} = obj.evaluate(pred_samples(:, i));
                end
                
            end
              
        end
        
        
        % evaluate an image star with specific values of predicates
        function image = evaluate(obj, pred_val)
            % @pred_val: valued vector of predicate variables
            
            % author: Dung Tran
            % date: 6/14/2019
            
            if isempty(obj.V)
                error('The imagestar is an empty set');
            end
            
            if size(pred_val, 2) ~= 1
                error('Invalid predicate vector');
            end
            
            if size(pred_val, 1) ~= obj.numPred
                error('Inconsistency between the size of the predicate vector and the number of predicates in the imagestar');
            end
            
            image(:, :, obj.numChannel) = zeros(obj.height, obj.width);
            
            for i=1:obj.numChannel
                
                image(:, :, i) = obj.V(:,:,i, 1);
                
                for j=2:obj.numPred + 1
                    
                    image(:, :, i) = image(:, :, i) + pred_val(j-1) * obj.V(:,:,i, j);
                
                end
                
            end
                      
        end
        
        
                
        % affineMap of an ImageStar is another imagestar
        % y = scale * x + offset;
        function image = affineMap(obj, scale, offset)
            % @scale: scale coefficient [1 x 1 x NumChannels] array
            % @offset: offset coefficient [1 x 1 x NumChannels] array
            % @image: a new ImageStar
            
            % author: Dung Tran
            % date: 1/1/2020
            
            
            if ~isempty(scale) && ~isscalar(scale) && size(scale, 3) ~= obj.numChannel
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
                         
            image = ImageStar(new_V, obj.C, obj.d, obj.pred_lb, obj.pred_ub);
                      
        end
        
        
        % transform to Star
        function S = toStar(obj)
            
            % author: Dung Tran
            % date: 7/19/2019
            
            nc = obj.numChannel;
            h = obj.height;
            w = obj.width;
            np = obj.numPred;
            
            N = h*w*nc; % total number of pixels in the input image         
            V1(:, np+1) = zeros(N, 1);
            for j=1:np+1
                V1(:,j) = reshape(obj.V(:,:,:, j), N, 1);
            end       
            if ~isempty(obj.im_lb) && ~isempty(obj.im_ub)
                state_lb = reshape(obj.im_lb, N, 1);
                state_ub = reshape(obj.im_ub, N, 1);
                S = Star(V1, obj.C, obj.d, obj.pred_lb, obj.pred_ub, state_lb, state_ub);
            else
                S = Star(V1, obj.C, obj.d, obj.pred_lb, obj.pred_ub);
            end
            
        end
        
        % checking if an ImageStar is an empty set
        function bool = isEmptySet(obj)
            % author: Dung Tran
            % date: 1/10/2020
            
            S = obj.toStar;
            bool = S.isEmptySet;
        end
        
        
        % contain, check if ImageStar contains an image
        function bool = contains(obj, image)
            % @image: input image
            % @bool: = 1 if the ImageStar contain the image
            %        = 0 if the ImageStar does not contain the image
            
            % author: Dung Tran
            % date: 1/8/2020
            
            n = size(image);
            if length(n) == 2 % one channel image
                if n(1) ~= obj.height || n(2) ~= obj.width || obj.numChannel ~= 1
                    error('Inconsistent dimenion between input image and the ImageStar');
                end
                y = reshape(image, [n(1)*n(2) 1]);
            elseif length(n) == 3
                if n(1) ~= obj.height || n(2) ~= obj.width || n(3) ~= obj.numChannel
                    error('Inconsistent dimenion between input image and the ImageStar');
                end
                y = reshape(image, [n(1)*n(2)*n(3) 1]);
            else
                error('Invalid input image');
            end
            
            S = obj.toStar;
            bool = S.contains(y);
            
        end
        
        
        
        % projection of imagestar on specific 2d plane
        function proj_star = project2D(obj, point1, point2)
            % @point1: index of first dimension
            % @point2: index of second dimension
            % @proj_star: is a projected star set
            
            
            % author: Dung Tran
            % date: 8/28/2019
            
            
            if (length(point1) ~= 3) || (length(point2) ~= 3)
                error('Invalid input point');
            end
            
            if ((point1(1) < 1) || (point1(1) > obj.height) || (point1(2) < 1) || (point1(2) > obj.width) || (point1(3) < 1) || (point1(3) > obj.numChannel))
                error('The first input point is invalid');
            end
            
            if ((point2(1) < 1) || (point2(1) > obj.height) || (point2(2) < 1) || (point2(2) > obj.width) || (point2(3) < 1) || (point2(3) > obj.numChannel))
                error('The first input point is invalid');
            end
            
            n = obj.numPred + 1; 
            
            new_V = zeros(2, n);
            for i=1:n
                new_V(1, i) = obj.V(point1(1), point1(2), point1(3), i);
                new_V(2, i) = obj.V(point2(1), point2(2), point2(3), i);              
            end     
            
            proj_star = Star(new_V, obj.C, obj.d, obj.pred_lb, obj.pred_ub);
            
        end
        
        
        % get ranges of a state at specific position
        function [xmin, xmax] = getRange(obj, vert_ind, horiz_ind, chan_ind)
            % @vert_ind: vectical index
            % @horiz_ind: horizontal index
            % @chan_ind : channel index
            % @xmin: min of x(vert_ind,horiz_ind, channel_ind)
            % @xmax: max of x(vert_ind,horiz_ind, channel_ind)
            
            
            % author: Dung Tran
            % date: 6/18/2019
                                  
            
            if isempty(obj.C) || isempty(obj.d)
                error('The imagestar is empty');
            end
            
            if vert_ind < 1 || vert_ind > obj.height
                error('Invalid veritical index');
            end
            
            if horiz_ind < 1 || horiz_ind > obj.width
                error('Invalid horizonal index');
            end
            
            if chan_ind < 1 || chan_ind > obj.numChannel
                error('Invalid channel index');
            end
            
               
            f = obj.V(vert_ind, horiz_ind, chan_ind, 2:obj.numPred + 1);
            [~, fval, exitflag, ~] = glpk(f, obj.C, obj.d, obj.pred_lb, obj.pred_ub);
            if exitflag == 5
                xmin = fval + obj.V(vert_ind, horiz_ind, chan_ind, 1);
            else
                error('Cannot find an optimal solution, exitflag = %d', exitflag);
            end          

            [~, fval, exitflag, ~] = glpk(-f, obj.C, obj.d, obj.pred_lb, obj.pred_ub);
            if exitflag == 5
                xmax = -fval + obj.V(vert_ind, horiz_ind, chan_ind, 1);
            else
                error('Cannot find an optimal solution exitflag = %d', exitflag);
            end
            
            obj.im_lb(vert_ind, horiz_ind, chan_ind) = xmin;
            obj.im_ub(vert_ind, horiz_ind, chan_ind) = xmax;
                   
        end
        
        
        % estimate range quickly using only predicate bound information
        function [xmin, xmax] = estimageRange(obj, h, w, c)
            % @h: height index
            % @w: width index
            % @c: channel index
            % @xmin: min of x[h, w, c]
            % @xmax: max of x[h, w, c]
            
            % author: Dung Tran
            % date: 7/22/2019
            
            if isempty(obj.C) || isempty(obj.d)
                error('The imagestar is empty');
            end
            
            if h < 1 || h > obj.height
                error('Invalid veritical index');
            end
            
            if w < 1 || w > obj.width
                error('Invalid horizonal index');
            end
            
            if c < 1 || c > obj.numChannel
                error('Invalid channel index');
            end
            
            f = obj.V(h, w, c, 1:obj.numPred + 1);
            xmin = f(1);
            xmax = f(1);
            
            for i=2:obj.numPred+1
                if f(i) >= 0
                    xmin = xmin + f(i) * obj.pred_lb(i-1);
                    xmax = xmax + f(i) * obj.pred_ub(i-1);
                else
                    xmin = xmin + f(i) * obj.pred_ub(i-1);
                    xmax = xmax + f(i) * obj.pred_lb(i-1);
                end
                
            end
            
            
        end
        
        
        % estimate ranges quickly using only predicate bound information
        function [image_lb, image_ub] = estimateRanges(obj)
            % @h: height index
            % @w: width index
            % @c: channel index
            % @image_lb: lower bound image
            % @image_ub: upper bound image
            
            % author: Dung Tran
            % date: 7/22/2019
            
            if isempty(obj.C) || isempty(obj.d)
                error('The imagestar is empty');
            end
            
            if isempty(obj.im_lb) || isempty(obj.im_ub)
                         
                image_lb = zeros(obj.height, obj.width, obj.numChannel);
                image_ub = zeros(obj.height, obj.width, obj.numChannel);

                for i=1:obj.height
                    for j=1:obj.width
                        for k=1:obj.numChannel                     
                            [image_lb(i, j, k), image_ub(i, j, k)] = obj.estimateRange(i,j,k);
                        end
                    end
                end

                obj.im_lb = image_lb;
                obj.im_ub = image_ub;
                
            else
                
                image_lb = obj.im_lb;
                image_ub = obj.im_ub;
                
            end
         
        end
        
        
        % get lowew bound and upper bound images of an imagestar
        function [image_lb, image_ub] = getRanges(obj)
            % @image_lb: lower bound image
            % @image_ub: upper bound image
            
            % author: Dung Tran
            % date: 6/20/2019
            
            image_lb = zeros(obj.height, obj.width, obj.numChannel);
            image_ub = zeros(obj.height, obj.width, obj.numChannel);
            
            for i=1:obj.height
                for j=1:obj.width
                    for k=1:obj.numChannel                     
                        [image_lb(i, j, k), image_ub(i, j, k)] = obj.getRange(i,j,k);
                    end
                end
            end
            
            obj.im_lb = image_lb;
            obj.im_ub = image_ub;
                   
        end
        
        
        % quickly estimate range
        function [xmin, xmax] = estimateRange(obj, vert_ind, horiz_ind, chan_ind)
            % @vert_ind: vectical index
            % @horiz_ind: horizontal index
            % @chan_ind : channel index
            % @xmin: min of x(vert_ind,horiz_ind, channel_ind)
            % @xmax: max of x(vert_ind,horiz_ind, channel_ind)
            
            
            % author: Dung Tran
            % date: 7/19/2019
                                  
            
            if isempty(obj.C) || isempty(obj.d)
                error('The imagestar is empty');
            end
            
            if vert_ind < 1 || vert_ind > obj.height
                error('Invalid veritical index');
            end
            
            if horiz_ind < 1 || horiz_ind > obj.width
                error('Invalid horizonal index');
            end
            
            if chan_ind < 1 || chan_ind > obj.numChannel
                error('Invalid channel index');
            end
            
               
            f = obj.V(vert_ind, horiz_ind, chan_ind, 1:obj.numPred + 1);
            xmin = f(1);
            xmax = f(1);
            
            for i=2:obj.numPred + 1
                if f(i) >= 0
                    xmin = xmin + f(i) * obj.pred_lb(i-1);
                    xmax = xmax + f(i) * obj.pred_ub(i-1);
                else
                    xmin = xmin + f(i) * obj.pred_ub(i-1);
                    xmax = xmax + f(i) * obj.pred_lb(i-1);
                end
            end
            
        end
        
        
        % update local ranges for Max Pooling operation
        function updateRanges(obj, points)
            % @points: local points = [x1 y1 c1; x2 y2 c2; ...]
            
            % author: Dung Tran
            % date: 6/25/2019
            

            n = size(points, 1);
            for i=1:n
                fprintf('\nUpdate range at point: [h = %d, w = %d, c = %d]', points(i, 1), points(i,2), points(i, 3));
                obj.getRange(points(i, 1), points(i,2), points(i, 3));
            end        
            
        end
        
        
        % get local bound for Max Pooling operation
        function [lb, ub] = get_localBound(obj, startpoint, PoolSize, channel_id)
            % @startpoint: startpoint of the local(partial) image
            %               startpoint = [x1 y1];
            % @PoolSize: = [height width] the height and width of max pooling layer
            % @channel_id: the index of the channel
            % @lb, ub: the lower bound and upper bound of all points in the
            % local region
            
            % author: Dung Tran
            % date: 6/25/2019
            
            points = obj.get_localPoints(startpoint, PoolSize);
            n = length(points);
            % get lower bound and upper bound image
            if isempty(obj.im_lb) || isempty(obj.im_ub)
                [image_lb, image_ub] = obj.getBox;
            else
                image_lb = obj.im_lb;
                image_ub = obj.im_ub;
            end
            
            lb = image_lb(points(1,1), points(1,2), channel_id);
            ub = image_ub(points(1,1), points(1,2), channel_id);
            
            for i=2:n
                if image_lb(points(i,1), points(i,2), channel_id) < lb
                    lb = image_lb(points(i,1), points(i,2), channel_id);
                end
                if image_ub(points(i,1), points(i,2), channel_id) > ub
                    ub = image_ub(points(i,1), points(i,2), channel_id);
                end
            end
            
            
        end
        
        % get all local points index for Max Pooling operation
        function points = get_localPoints(obj, startpoint, PoolSize)
            % @startpoint: startpoint of the local(partial) image
            %               startpoint = [x1 y1];
            % @PoolSize: = [height width] the height and width of max pooling layer
            % @points: all indexes of all points for a single max
            % pooling operation (including the startpoint)
            
            % author: Dung Tran
            % date: 6/25/2019
            
            x0 = startpoint(1); % vertical index of the startpoint
            y0 = startpoint(2); % horizontal index of the startpoint
            h  = PoolSize(1);   % height of the MaxPooling layer
            w  = PoolSize(2);   % width of the MaxPooling layer
            

            if x0 < 1 || y0 < 1 || x0 + h - 1 > obj.height || y0 + w - 1 > obj.width
                error('Invalid startpoint or PoolSize');
            end
            points = zeros(h*w, 2);
            for i=1:h
                if i==1
                    x1 = x0;
                else
                    x1 = x1 + 1;
                end
                
                for j=1:w
                    if j==1
                        y1 = y0;
                    else
                        y1 = y1 + 1;
                    end
                    points((i-1)*w + j, :) = [x1 y1];
                end
            end
                    
        end
            
                    
        % get local max index, this medthod tries to find the maximum point
        % of a local image, used in over-approximate reachability analysis
        % of maxpooling operation
        function max_id = get_localMax_index(obj, startpoint, PoolSize, channel_id)
            % @startpoint: startpoint of the local(partial) image
            %               startpoint = [x1 y1];
            % @PoolSize: = [height width] the height and width of max pooling layer
            % @channel_id: the channel index
            % @max_id: = []: we don't know which one has maximum value,
            % i.e., the maximum values may be the intersection between of
            % several pixel valutes.
            %           = [xi yi]: the point that has maximum value
            
            % author: Dung Tran
            % date: 6/24/2019
            
            points = obj.get_localPoints(startpoint, PoolSize);          
            % get lower bound and upper bound image
            if isempty(obj.im_lb) || isempty(obj.im_ub)
                obj.estimateRanges;
            end
            
            h  = PoolSize(1);   % height of the MaxPooling layer
            w  = PoolSize(2);   % width of the MaxPooling layer
            n = h*w;
            
            lb = zeros(1, n);
            ub = zeros(1,n);
            
            for i=1:n
                point_i = points(i, :);
                lb(:, i) = obj.im_lb(point_i(1), point_i(2), channel_id);
                ub(:, i) = obj.im_ub(point_i(1), point_i(2), channel_id);   
            end
            
            [max_lb_val, max_lb_idx] = max(lb, [], 2);
            
            a = find(ub - max_lb_val > 0);
            a1 = find(ub - max_lb_val >= 0);
            a(a==max_lb_idx) = [];
                      
            if isempty(a)              
                max_id = points(max_lb_idx, :);
            else
                candidates = a1;
                % update local ranges               
                m = length(candidates);
                new_points = zeros(m, 3);
                new_points1 = zeros(m, 2);
                for i=1:m
                    p = points(candidates(i), :);
                    new_points = [p channel_id];
                    new_points1(i, :) = p;
                end
                obj.updateRanges(new_points);
                                
                lb = zeros(1, m);
                ub = zeros(1, m);

                for i=1:m
                    point_i = points(candidates(i), :);
                    lb(:, i) = obj.im_lb(point_i(1), point_i(2), channel_id);
                    ub(:, i) = obj.im_ub(point_i(1), point_i(2), channel_id);   
                end

                [max_lb_val, max_lb_idx] = max(lb, [], 2);

                a = find(ub - max_lb_val > 0);
                a(a==max_lb_idx) = [];
                                               
                if isempty(a) == 1
                    max_id = new_points1(max_lb_idx, :);
                else
                    candidates1 = find(ub - max_lb_val >= 0);
                    max_id = new_points1(max_lb_idx, :);
                    candidates1(candidates1 == max_lb_idx) = [];                    
                    m = length(candidates1);
                    max_id1 = max_id;
                    for j=1:m
                        p1 = new_points1(candidates1(j), :);         
                        if obj.is_p1_larger_p2([p1(1) p1(2) channel_id], [max_id(1) max_id(2) channel_id])
                            max_id1 = [max_id1; p1];
                        end
                    end                
                    max_id = max_id1; 
                    fprintf('\nThe local image has %d max candidates.', size(max_id,1));
                    
                end
                
                n = size(max_id, 1);
                channels = channel_id*ones(n,1);
                max_id = [max_id channels];

            end
         
         
        end
        
        
        
        % compare between two specific points in the image
        % check if p1 > p2 is feasible or not
        % useful for max pooling operation
        function b = is_p1_larger_p2(obj, p1, p2)
            % @p1: the first point = [h1, w1, c1]
            % @p2: the second point = [h2, w2, c2]
            % h: height, w: width, c: channel index
            
            % @b = 1 -> p1 > p2 is feasible
            %    = 0 -> p1 > p2 is not feasible
            
            % author: Dung Tran
            % date: 7/23/2019
            
            % a*(V2 - V1) <= c1 - c2
            
            C1 = zeros(1, obj.numPred);
            for i=2:obj.numPred+1
                C1(i-1) = obj.V(p2(1),p2(2),p2(3), i) - obj.V(p1(1), p1(2), p1(3), i);
            end
            
            d1 = obj.V(p1(1), p1(2), p1(3), 1) - obj.V(p2(1), p2(2), p2(3)); 
            
            new_C = [obj.C; C1];
            new_d = [obj.d; d1];
           
            f = zeros(1, obj.numPred);

                [~,~,status,~] = glpk(f, new_C, new_d, obj.pred_lb, obj.pred_ub);

                if status == 5 % feasible solution exist
                    b = 1;
                else
                    b = 0;
                end
            
        end
               
              
         
    end
    
    
    
    methods(Static)
        
        
        % check if a pixel value is the maximum value compared with others
        % this is core step for exactly performing maxpooling operation on an
        % imagestar set
        function [new_C, new_d] = isMax(maxMap, ori_image, center, others)
            % @maxMap: the current maxMap ImageStar
            % @ori_image: the original ImageStar to compute the maxMap 
            % @center: is the center pixel position we want to check
            %          center = [x1 y1 c1]
            % @others: is the other pixel position we want to compare with
            % the center one
            %          others = [x2 y2 c2; x3 y3 c3]
            % @out_image: = imagestar object = in_image with with some updates in the predicate
            %         constraints
            
            % author: Dung Tran
            % date: 6/20/2019
            % update: 7/25/2019
                       
            if maxMap.numPred ~= ori_image.numPred
                error('Inconsistency between number of predicates in the current maxMap and the original image');
            end
            
            n = size(others, 1);
           
            % the center may be the max point with some extra
            % constraints on the predicate variables
            new_C = zeros(n, maxMap.numPred);
            new_d = zeros(n, 1);
            
            
            for i=1:n                
                % add new constraint
                % compare point (i,j) with point (i1, j1)
                % p[i,j] = c[i,j] + \Sigma (V^k[i,j]*a_k), k=1:numPred
                % p[i1,j1] = c[i1, j1] + \Sigma ((V^k[i1,j1]*a_k), k=1:numPred)

                % p[i,j] >= p[i1, j1] <=> 
                % <=> \Sigma (V^k[i1, j1] - V^k[i,j])*a[k] <= c[i,j] - c[i1,j1]

                new_d(i) = ori_image.V(center(1),center(2), center(3), 1) - ori_image.V(others(i,1), others(i,2), others(i,3), 1);
                for j=1:maxMap.numPred                    
                    new_C(i,j) = -ori_image.V(center(1),center(2), center(3), j+1) + ori_image.V(others(i,1), others(i,2), others(i,3), j+1);
                end           
            end
            
            
            C1 = [maxMap.C; new_C];
            d1 = [maxMap.d; new_d];

            % remove redundant constraints
            E = [C1 d1];
            E = unique(E, 'rows');

            C1 = E(:, 1:ori_image.numPred);
            d1 = E(:, ori_image.numPred + 1);

            f = zeros(1, ori_image.numPred);

            [~,~,status,~] = glpk(f, C1, d1, ori_image.pred_lb, ori_image.pred_ub);

            if status == 5 % feasible solution exist
                new_C = C1;
                new_d = d1;
            else
                new_C = [];
                new_d = [];
            end

        end
        
        % step split of an image star
        % a single in_image can be splitted into several images in the
        % exact max pooling operation
        function images = stepSplit(in_image, ori_image, pos, split_index)
            % @in_image: the current maxMap ImageStar
            % @ori_image: the original ImageStar to compute the maxMap 
            % @pos: local position of the maxMap where splits may occur
            % @split_index: indexes of local pixels where splits occur
            
            % author: Dung Tran
            % date: 7/25/2019
            
            
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
                [new_C, new_d] = ImageStar.isMax(in_image, ori_image, center, others);                
                if ~isempty(new_C) && ~isempty(new_d)                    
                    V = in_image.V;
                    V(pos(1), pos(2), pos(3), :) = ori_image.V(center(1), center(2), center(3), :);
                    im = ImageStar(V, new_C, new_d, in_image.pred_lb, in_image.pred_ub, in_image.im_lb, in_image.im_ub);
                    images = [images im];
                end
            end
           
        end
        
        
        % step split for multiple image stars
        % a single in_image can be splitted into several images in the
        % exact max pooling operation
        function images = stepSplitMultipleInputs(in_images, ori_image, pos, split_index, option)
            % @in_image: the current maxMap ImageStar
            % @ori_image: the original ImageStar to compute the maxMap 
            % @pos: local position of the maxMap where splits may occur
            % @split_index: indexes of local pixels where splits occur
            % @option: = [] or 'parallel'
            
            % author: Dung Tran
            % date: 7/25/2019
            
            
            n = length(in_images);
            images = [];
            if strcmp(option, 'parallel')
                parfor i=1:n
                    images = [images ImageStar.stepSplit(in_images(i), ori_image, pos, split_index)];
                end
            elseif isempty(option) || strcmp(option, 'single')
                for i=1:n
                    images = [images ImageStar.stepSplit(in_images(i), ori_image, pos, split_index)];
                end
            else 
                error('Unknown computation option');
            end       
            
        end
        
        
        
        
        % reshape an ImageStar
        function new_IS = reshape(in_IS, new_shape)
            % @in_IS: input ImageStar
            % @new_shape: new shape
            % @new_IS: new ImageStar

            % author: Dung Tran
            % date: 8/21/2019

            n = size(new_shape);

            if n(2) ~= 3 || n(1) ~= 1
                error('new shape should be 1x3 array');
            end

            if new_shape(1) * new_shape(2) * new_shape(3) ~= in_IS.height * in_IS.width * in_IS.numChannel
                error('new shape is inconsistent with the current shape');
            end
            
            new_V = reshape(in_IS.V, [new_shape in_IS.numPred + 1]);           
            new_IS = ImageStar(new_V, in_IS.C, in_IS.d, in_IS.pred_lb, in_IS.pred_ub, in_IS.im_lb, in_IS.im_ub);
     
        end
        
        % add new constraint to predicate variables of an ImageStar
        % used for finding counter examples
        % we will add new constraint: p2 >= p1
        function [new_C, new_d] = addConstraint(in_IS, p1, p2)
            % @in_IS: input ImageStar
            % @p1: first point position
            % @p2: second point position
            % @new_C: a new predicate constraint matrix C 
            % @new_d: new predicate constraint vector
            
            % author: Dung Tran
            % date: 8/21/2019
            
            
            if ~isa(in_IS, 'ImageStar')
                error('Input set is not an ImageStar');
            end
            
            new_d = in_IS.V(p2(1), p2(2), p2(3), 1) - in_IS.V(p1(1), p1(2), p1(3), 1);
            new_C = in_IS.V(p1(1), p1(2), p1(3), 2:in_IS.numPred + 1) - in_IS.V(p2(1), p2(2), p2(3), 2:in_IS.numPred + 1);
            
            new_C = [in_IS.C; new_C];
            new_d = [in_IS.d; new_d];
            
        end



        
        
    end
    
    
    
    
    
    
   
        
       
end

