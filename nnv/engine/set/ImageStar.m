classdef ImageStar
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
        V = []; % a cell (size = numChannel) of 2-d basis matrices 
        C = []; % a constraints matrix of the predicate
        d = []; % a constraints vector of the predicate
        numPred = 0; % number of predicate variables
        pred_lb = []; % lower bound vector of the predicate
        pred_ub = []; % upper bound vector of the predicate
        im_lb = []; % lower bound image of the ImageStar
        im_ub = []; % upper bound image of the ImageStar
        
        % used to store max points in max pooling analysis
        % a single max_point = [height_position; width_position; channel_index] 
        max_points = [];        

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
                    
                    obj.numPred = 0;
                    C = [];
                    obj.d = [];
                    obj.pred_lb = [];
                    obj.pred_ub = [];
                    
                    % get number of predicate variables
                    for k=1:obj.numChannel
                        for i=1:obj.height
                            for j=1:obj.width
                                if obj.LB(i,j,k) < obj.UB(i, j, k)
                                    obj.numPred = obj.numPred + 1;
                                end
                            end
                        end
                    end
                    
                    V1(:, :, obj.numPred, obj.numChannel) = zeros(obj.height, obj.width);
                    
                    pred_count = 0;
                    for i=1:obj.numChannel
                        c = reshape(obj.IM(:,:,i)', [obj.height * obj.width,1]);
                        lb = reshape(obj.LB(:,:,i)', [obj.height * obj.width,1]);
                        ub = reshape(obj.UB(:,:,i)', [obj.height * obj.width,1]);
                        lb = lb + c;
                        ub = ub + c;
                        B = Box(lb, ub);
                        X = B.toStar;                  
                                                                       
                        for j=1:X.nVar + 1    
                            A = reshape(X.V(:,j), [obj.height, obj.width]);
                            if j==1
                                V1(:,:,1, i) = A';
                            else
                                V1(:, :, j+pred_count, i) = A';
                            end
                        end 
                        pred_count = pred_count + X.nVar;

                        obj.C = blkdiag(obj.C, X.C);
                        obj.d = [obj.d; X.d];
                        obj.pred_lb = [obj.pred_lb; X.predicate_lb];
                        obj.pred_ub = [obj.pred_ub; X.predicate_ub];                                                     
                    end
                    
                    obj.V = V1;
                                         
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
                        obj.numChannel = 1;
                    elseif length(n) == 4
                        obj.numChannel = n(4);
                    else
                        error('Invalid basis matrix');
                    end
                                            
                    obj.V = V1;
                    obj.height = n(1);
                    obj.width = n(2);
                    
                    
                case 8 % 
                    
                    V1 = varargin{1};  % basis matrices
                    C1 = varargin{2};  % predicate constraint matrix 
                    d1 = varargin{3};  % predicate constraint vector
                    lb1 = varargin{4}; % predicate lower bound
                    ub1 = varargin{5}; % predicate upper bound
                    im_lb1 = varargin{6}; % lower bound image
                    im_ub1 = varargin{7}; % upper bound image
                    mp = varargin{8}; % maxpoints
                    
                    
                                        
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
                        obj.numChannel = 1;
                    elseif length(n) == 4
                        obj.numChannel = n(4);
                    else
                        error('Invalid basis matrix');
                    end
                                            
                    obj.V = V1;
                    obj.height = n(1);
                    obj.width = n(2);
                    
                    if size(im_lb1,1) ~= obj.height || size(im_lb1, 2) ~= obj.width
                        error('Inconsistent dimension between lower bound image and the constructed imagestar');
                    else
                        obj.im_lb = im_lb1;
                    end
                    
                    if size(im_ub1,1) ~= obj.height || size(im_ub1, 2) ~= obj.width
                        error('Inconsistent dimension between upper bound image and the constructed imagestar');
                    else
                        obj.im_ub = im_ub1;
                    end
                    
                    if ~isempty(mp) && size(mp, 1) ~= 3
                        error('Invalid max_points matrix');
                    else
                        obj.max_points = mp;
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
                    
                    error('Invalid number of input arguments, (should be from 0, 3, 5 , 8)');
                    
            end
                        
        end
        
        
        % add max points
        function obj = add_maxPoint(obj, height_pos, width_pos, channel_ind)
            % @height_pos: height position of the max point
            % @width_pos: width position of the max point
            % @channel_ind: index of the channel that tha max point belongs
            % to 
            
            % author: Dung Tran
            % date: 6/20/2019
            
            
            if height_pos < 1 || height_pos > obj.height
                error('Invalid height position');
            end
            
            if width_pos < 1 || width_pos > obj.width
                error('Invalid width position');
            end
            
            if channel_ind < 1 || channel_ind > obj.numChannel
                error('Invalid channel index');
            end
            
            new_point = [height_pos; width_pos; channel_ind];
            obj.max_points = [obj.max_points new_point];
            
        end
        
        
        % randomly show images sampled from image star
        function show(obj, N)
            % @N: number of sampled images 

            % author: Dung Tran
            % date: 6/14/2019
            
            
            if isempty(obj.V)
                error('The imagestar is an empty set');
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
                
                image(:, :, i) = obj.V(:,:,1,i);
                
                for j=2:obj.numPred + 1
                    
                    image(:, :, i) = image(:, :, i) + pred_val(j-1) * obj.V(:,:,j,i);
                
                end
                
            end
                      
        end
        
        
        % get ranges of a state at specific position
        function [xmin, xmax] = getRange(varargin)
            % @vert_ind: vectical index
            % @horiz_ind: horizontal index
            % @xmin: min of x(vert_ind,horiz_ind)
            % @xmax: max of x(vert_ind,horiz_ind)
            % @option: = 'parallel' use parallel computing
            
            % author: Dung Tran
            % date: 6/18/2019
            
            switch nargin
                case 4
                    obj = varargin{1};
                    vert_ind = varargin{2};
                    horiz_ind = varargin{3};
                    option = varargin{4};
                    
                    if ~strcmp(option, 'parallel')
                        error('Unknown option for computation');
                    end
                
                case 3
                    obj = varargin{1};
                    vert_ind = varargin{2};
                    horiz_ind = varargin{3};
                    option = [];
                    
                otherwise
                    error('Invalid number of arguments, should be 2, or 3');
            end
            
            
            if isempty(obj.C) || isempty(obj.d)
                error('The imagestar is empty');
            end
            
            if vert_ind < 1 || vert_ind > obj.height
                error('Invalid veritical index');
            end
            
            if horiz_ind < 1 || horiz_ind > obj.width
                error('Invalid horizonal index');
            end
            
            xmin = zeros(obj.numChannel,1);
            xmax = zeros(obj.numChannel,1);

            if strcmp(option, 'parallel')
                
                parfor i=1:obj.numChannel
                
                    f = obj.V(vert_ind, horiz_ind, 2:obj.numPred + 1, i);
                    [~, fval, exitflag, ~] = glpk(f, obj.C, obj.d);
                    if exitflag == 5
                        xmin(i) = fval + obj.V(vert_ind, horiz_ind, 1, i);
                    else
                        error('Cannot find an optimal solution, exitflag = %d', exitflag);
                    end          

                    [~, fval, exitflag, ~] = glpk(-f, obj.C, obj.d);
                    if exitflag > 0
                        xmax(i) = -fval + obj.V(vert_ind, horiz_ind, 1, i);
                    else
                        error('Cannot find an optimal solution exitflag = %d', exitflag);
                    end
                             
                end
                
            else
                
                for i=1:obj.numChannel
                    f = obj.V(vert_ind, horiz_ind, 2:obj.numPred + 1, i);
                    [~, fval, exitflag, ~] = glpk(f, obj.C, obj.d);
                    if exitflag == 5
                        xmin(i) = fval + obj.V(vert_ind, horiz_ind, 1, i);
                    else
                        error('Cannot find an optimal solution, exitflag = %d', exitflag);
                    end          

                    [~, fval, exitflag, ~] = glpk(-f, obj.C, obj.d);
                    if exitflag == 5
                        xmax(i) = -fval + obj.V(vert_ind, horiz_ind, 1, i);
                    else
                        error('Cannot find an optimal solution exitflag = %d', exitflag);
                    end

                end
                
            end    
            
        end
        
        
        % get lowew bound and upper bound images of an imagestar
        function [image_lb, image_ub] = getBox(obj)
            % @image_lb: lower bound image
            % @image_ub: upper bound image
            
            % author: Dung Tran
            % date: 6/20/2019
            
            image_lb = zeros(obj.height, obj.width, obj.numChannel);
            image_ub = zeros(obj.height, obj.width, obj.numChannel);
            
            for i=1:obj.height
                for j=1:obj.width
                    [xmin, xmax] = obj.getRange(i,j);
                    %[xmin, xmax] = obj.getRange(i,j, 'parallel');
                    image_lb(i, j, :) = xmin';
                    image_ub(i,j,:) = xmax';
                end
            end
            
            
        end
             
                
        % check if a pixel value is the maximum value compared with others
        % this is core step for exactly performing maxpooling operation on an
        % imagestar set
        function image = isMax(obj, center, others, channel_ind)
            % @center: is the center pixel position we want to check
            %          center = [i; j]
            % @others: is the other pixel position we want to compare with
            % the center one
            %          others = [i1 i2 ... in, j1 j2 ... jn]
            %          others(:, i) is the position of i^th pixel
            % @channel_ind: the index of current channel
            % @image: = [], the pixel value cannot be the maximum one
            %         = imagestar object, the pixel value can be the
            %         maximum one with some updates in the predicate
            %         constraints
            
            % author: Dung Tran
            % date: 6/20/2019
            
            
            if length(center) ~= 2 || center(1) > obj.height || center(1) < 1 || center(2) > obj.width || center(2) < 1
                error('Invalid center point');
            end
            
            if channel_ind > obj.numChannel || channel_ind < 1
                error('Invalid channel index');
            end
            
            if size(others, 1) ~= 2
                error('Invalid position matrix');
            end
            
            if isempty(obj.im_lb) || isempty(obj.im_ub)
                [obj.im_lb, obj.im_ub] = obj.getBox;
            end
            n = size(others, 2);
            min_center = obj.im_lb(center(1), center(2), channel_ind);
            max_center = obj.im_ub(center(1), center(2), channel_ind);
            max_count = 0;
            min_count = 0;
            for i=1:n
                min_point_i = obj.im_lb(others(1,i), others(2,i), channel_ind);
                max_point_i = obj.im_lb(others(1, i), others(2, i), channel_ind);
                if min_center >= max_point_i
                    max_count = max_count + 1;
                end
                if max_center <= min_point_i
                    min_count = min_count + 1;
                end
            end
            
            if max_count == n % the center is the max point by only checking the bound
                image = ImageStar(obj.V, obj.C, obj.d, obj.pred_lb, obj.pred_ub, obj.im_lb, obj.im_ub, obj.max_points);
                image = image.add_maxPoint(center(1), center(2), channel_ind);
            elseif min_count == n % the center is not the max point by only checking the bound
                image = [];
            else
                % the center may be the max point with some extra
                % constraints on the predicate variables
                new_C = zeros(n, obj.numPred);
                new_d = zeros(n, 1);

                for i=1:n                
                    % add new constraint
                    % compare point (i,j) with point (i1, j1)
                    % p[i,j] = c[i,j] + \Sigma (V^k[i,j]*a_k), k=1:numPred
                    % p[i1,j1] = c[i1, j1] + \Sigma ((V^k[i1,j1]*a_k), k=1:numPred)

                    % p[i,j] >= p[i1, j1] <=> 
                    % <=> \Sigma (V^k[i1, j1] - V^k[i,j])*a[k] <= c[i,j] - c[i1,j1]

                    new_d(i) = obj.V(center(1),center(2), 1, channel_ind) - obj.V(others(1,i), others(2,i), 1, channel_ind);
                    for j=1:obj.numPred                    
                        new_C(i,j) = obj.V(center(1),center(2), j+1, channel_ind) - obj.V(others(1,i), others(2,i), j+1, channel_ind);
                    end           
                end

                C1 = [obj.C; new_C];
                d1 = [obj.d; new_d];

                % remove redundant constraints
                E = [C1 d1];
                E = unique(E, 'rows');

                C1 = E(:, 1:obj.numPred);
                d1 = E(:, obj.numPred + 1);

                f = zeros(1, obj.numPred);

                [~,~,status,~] = glpk(f, C1, d1);

                if status == 5 % feasible solution exist
                    image = ImageStar(obj.V, C1, d1, obj.pred_lb, obj.pred_ub, obj.im_lb, obj.im_ub, obj.max_points);
                    image = image.add_maxPoint(center(1), center(2), channel_ind);
                else
                    image = [];
                end
                
            end
            
        end
        
        
        
    end
    
    
    
   
        
       
end

