classdef VolumeStar < handle
    % Class for representing set of 3D/volume images using Star set
    % 
    % A volume image can be attacked by bounded noise. An attacked volume 
    % image can be represented using an VolumeStar Set
    % 
    % Diego Manzanas Lopez: October 12, 2023
    % 
    %=================================================================%
    %   A volume color image is represented by 4-dimensional array.
    % 
    %   A volume image is a h x w x dp x c matrix, 
    %    - h: height of volume image
    %    - w: width of volume image
    %    - dp: depth of volume image
    %    - c: channels of volume image
    %   
    % *** A gray-scale volume image has only one channel.
    %
    %=================================================================%
    
    
    properties
        numChannel = 0; % number of channels
        height = 0;     % height of volume image
        width = 0;      % width of volume image
        depth = 0;      % depth of volume image
        
        % These 3 variables are useful for initializing the input set to NN
        IM = []; % center volume image (high-dimensional array)
        LB = []; % lower bound of attack (high-dimensional array)
        UB = []; % upper bound of attack (high-dimensional array)
        
        % These variables are used to represent the set, used for
        % reachability analysis
        V = []; % a cell (size = numPred)
        C = []; % a constraints matrix of the predicate
        d = []; % a constraints vector of the predicate
        numPred = 0; % number of predicate variables
        pred_lb = []; % lower bound vector of the predicate
        pred_ub = []; % upper bound vector of the predicate
        vol_lb = []; % lower bound volume of the VolumeStar
        vol_ub = []; % upper bound bolume of the VolumeStar
        
        MaxIdxs = cell(1,1); % used for unmaxpooling operation in Segmentation network
        InputSizes = cell(1,1); % used for unmaxpooling operation in Segmentation network
        
    end
    
    methods % Constructor and sampling methods

        % constructor using variables or bounds of a VolumeStar
        function obj = VolumeStar(varargin)
            % @nargin = 3: IM = varargin{1}, LB = varagin{2}, UB =
            % varargin{3}
            %         = 2: Stars = varargin{1}, volumeSize = varagin{2}
            %         = otherwise: IM = [], LB = [], UB = [], Stars = []
            % @IM: center volume image (high-dimensional array)
            % @LB: lower bound of attack (high-dimensional array)
            % @UB: upper bound of attack (high-dimensional array)
            % @Stars: 1D representation of a VolumeStar (flattened VolumeStar)
            
            switch nargin
                
                case 3 % input center volume image and lower and upper bound matrices (box-representation)
                    
                    IM1 = varargin{1};
                    LB1 = varargin{2};
                    UB1 = varargin{3};
                    n = size(IM1); % n(1), n(2), n(3) are height, width and depth of volume image
                                   % n(4) is number of channels
                    l = size(LB1);
                    u = size(UB1);
                    
                    if n(1) ~= l(1) || n(1) ~= u(1) || n(2) ~= l(2) || n(2) ~= u(2) || n(3) ~= l(3) || n(3) ~= u(3) 
                        error('Inconsistency between center volume image and attack bound matrices');
                    end

                    if length(n) ~= length(l) || length(n) ~= length(u)
                        error('Inconsistency between center volume image and attack bound matrices');
                    end

                    if length(n) == 3 && length(l) == 3 && length(u) == 3
                        obj.numChannel = 1;
                        obj.IM = IM1;
                        obj.LB = LB1;
                        obj.UB = UB1;
                        obj.height = n(1);
                        obj.width = n(2);
                        obj.depth = n(3);
                    elseif length(n) == 4 && length(l) == 4 && length(u) == 4
                        if n(4) == l(4) && n(4) == u(4)
                            obj.numChannel = n(4);
                            obj.IM = IM1; 
                            obj.LB = LB1;
                            obj.UB = UB1;
                            obj.height = n(1);
                            obj.width = n(2);
                            obj.depth = n(3);
                        else
                            error('Inconsistent number of channels between the center volume image and the bound matrices');
                        end
                    else
                        error('Inconsistent number of channels between the center volume image and the bound matrices');
                    end
                    
                    obj.vol_lb = IM1 + LB1; % lower bound volume
                    obj.vol_ub = IM1 + UB1; % upper bound volume

                    % converting box VolumeStar to an array of 2D Stars
                    n = size(obj.vol_lb);
                    I = Star(reshape(obj.vol_lb, [prod(n), 1]), reshape(obj.vol_ub, [prod(n), 1]));
                    if length(n) == 4
                        obj.V = reshape(I.V,[n(1), n(2), n(3), n(4), I.nVar + 1]);
                    else
                        obj.V = reshape(I.V,[n(1), n(2), n(3),  1,   I.nVar + 1]);
                    end
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
                    
                    if length(n) == 4
                        obj.V = V1;
                        obj.height = n(1);
                        obj.width  = n(2);
                        obj.depth  = n(3);
                        obj.numChannel = n(4);
                    elseif length(n) == 5
                        if n(5) ~= obj.numPred + 1
                            error('Inconsistency between the basis matrix and the number of predicate variables');
                        else
                            obj.numChannel = n(4);
                            obj.V = V1;
                            obj.height = n(1);
                            obj.width  = n(2);
                            obj.depth  = n(3);
                        end
                    elseif length(n) == 3
                            obj.numChannel = 1;
                            obj.V = V1;
                            obj.height = n(1);
                            obj.width  = n(2);
                            obj.depth  = n(3);
                    else
                        error('Invalid basis matrix');
                    end
                                            
                case 7 % 
                    
                    V1 = varargin{1};  % basis matrices
                    C1 = varargin{2};  % predicate constraint matrix 
                    d1 = varargin{3};  % predicate constraint vector
                    lb1 = varargin{4}; % predicate lower bound
                    ub1 = varargin{5}; % predicate upper bound
                    vol_lb1 = varargin{6}; % lower bound volume
                    vol_ub1 = varargin{7}; % upper bound volume
                                                          
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
                    
                    if length(n) == 4
                        obj.V = V1;
                        obj.height = n(1);
                        obj.width  = n(2);
                        obj.depth  = n(3);
                        obj.numChannel = n(4);
                    elseif length(n) == 5
                        if n(5) ~= obj.numPred + 1
                            error('Inconsistency between the basis matrix and the number of predicate variables');
                        else
                            obj.numChannel = n(4);
                            obj.V = V1;
                            obj.height = n(1);
                            obj.width  = n(2);
                            obj.depth  = n(3);
                        end
                    elseif length(n) == 3
                            obj.numChannel = 1;
                            obj.V = V1;
                            obj.height = n(1);
                            obj.width  = n(2);
                            obj.depth  = n(3);
                    else
                        error('Invalid basis matrix');
                    end
                    
                    if ~isempty(vol_lb1) && (size(vol_lb1,1) ~= obj.height || size(vol_lb1, 2) ~= obj.width) || size(vol_lb1, 3) ~= obj.depth
                        error('Inconsistent dimension between lower bound volume and the constructed VolumeStar');
                    else
                        obj.vol_lb = vol_lb1;
                    end
                    
                    if ~isempty(vol_ub1) && (size(vol_ub1,1) ~= obj.height || size(vol_ub1, 2) ~= obj.width) || size(vol_ub1, 3) ~= obj.depth
                        error('Inconsistent dimension between upper bound volume and the constructed VolumeStar');
                    else
                        obj.vol_ub = vol_ub1;
                    end
                    
                case 2
                    % input lower bound and upper bound volume
                    lb_im = varargin{1};
                    ub_im = varargin{2};
                    if isempty(lb_im)
                        error('Invalid lower bound volume');
                    end
                    if isempty(ub_im)
                        error('Invalid upper bound volume');
                    end
                    
                    n = size(lb_im);
                    m = size(ub_im); 
                    
                    if length(n) ~= length(m)
                        error('Inconsistency between lower bound volume and upper bound volume.');
                    end
                    
                    if length(n) == 4
                        if n(1) ~= m(1) || n(2) ~= m(2) || n(3) ~= m(3) || n(4) ~= m(4)
                            error('Inconsistency between lower bound volume and upper bound volume.');
                        end
                        lb = reshape(lb_im, [prod(n) 1]);
                        ub = reshape(ub_im, [prod(n) 1]);
                        S = Star(lb, ub);
                        obj = S.toVolumeStar(n(1), n(2), n(3), n(4));
                        obj.vol_lb = lb_im;
                        obj.vol_ub = ub_im;
                    elseif length(n) == 3
                        if n(1) ~= m(1) || n(2) ~= m(2) || n(3) ~= m(3)
                            error('Inconsistency between lower bound volume and upper bound volume');
                        end
                        obj.numChannel = 1;
                        lb = reshape(lb_im, [prod(n) 1]);
                        ub = reshape(ub_im, [prod(n) 1]);
                        S = Star(lb, ub);
                        obj = S.toVolumeStar(n(1), n(2), n(3), 1);
                        obj.vol_lb = lb_im;
                        obj.vol_ub = ub_im;
                    end
                                         
                case 0 % create an empty VolumeStar
                                      
                otherwise
                    error('Invalid number of input arguments, (should be from 0, 3, 5 , or 7)');
            end
                        
        end
                
        % randomly generate a set of volumes from a VolumeStar set
        function volumes = sample(obj, N)
            % @N: number of volume images
            
            if isempty(obj.V)
                error('The VolumeStar is an empty set');
            end
            
            if isempty(obj.C) || isempty(obj.d)
                volumes = obj.IM;
            else
                V1 = [cast(zeros(obj.numPred, 1),'like', obj.C)  eye(obj.numPred)];
                S = Star(V1, obj.C, obj.d, obj.pred_lb, obj.pred_ub);
                pred_samples = S.sample(N); 
                
                M = size(pred_samples,2);
                volumes = cell(1, M);
                for i=1:M
                    volumes{i} = obj.evaluate(pred_samples(:, i));
                end
            end
              
        end
        
        % evaluate a VolumeStar with specific values of predicates
        function volume = evaluate(obj, pred_val)
            % @pred_val: valued vector of predicate variables
            
            if isempty(obj.V)
                error('The VolumeStar is an empty set');
            end
            
            if size(pred_val, 2) ~= 1
                error('Invalid predicate vector');
            end
            
            if size(pred_val, 1) ~= obj.numPred
                error('Inconsistency between the size of the predicate vector and the number of predicates in the VolumeStar');
            end
            
            volume(:, :, :, obj.numChannel) = cast(zeros(obj.height, obj.width, obj.depth), 'like', obj.V);
            
            for i=1:obj.numChannel
                volume(:, :, :, i) = obj.V(:, :, :, i, 1);
                for j=2:obj.numPred + 1
                    volume(:, :, i) = volume(:, :, :, i) + pred_val(j-1) * obj.V(:, :, :, i, j);
                end
            end
                      
        end
    
    end


    methods % operations / transformation functions

        % affineMap of a VolumeStar is another VolumeStar (y = scale * x + offset)
        function volume = affineMap(obj, scale, offset)
            % @scale: scale coefficient [1 x 1 x NumChannels] array
            % @offset: offset coefficient [1 x 1 x NumChannels] array
            % @volume: a new VolumeStar
            
            if ~isempty(scale) && ~isscalar(scale) && size(scale, 4) ~= obj.numChannel
                error('Inconsistent number of channels between scale array and the VolumeStar');
            end
            
            if ~isempty(scale) && (isvector(scale) || isscalar(scale))
                new_V = scale.*obj.V;
            elseif ~isempty(scale) && ismatrix(scale)
                new_V = pagemtimes(scale,obj.V);
            else
                new_V = obj.V;
            end
            
            if ~isempty(offset)
                new_V(:,:,:,:,1) = new_V(:,:,:,:,1) + offset;
            end
            volume = VolumeStar(new_V, obj.C, obj.d, obj.pred_lb, obj.pred_ub);
                      
        end
        
        % Minkowski Sum of two VolumeStars is another VolumeStar (y = x1 + x2)
        function volume = MinkowskiSum(obj, I)
            
            S1 = obj.toStar;
            S2 = I.toStar;
            S = S1.MinkowskiSum(S2);
            
            volume = S.toVolumeStar(I.height, I.width, I.depth, I.numChannel); 

        end
        
        % concatenation of two VolumeStars is another VolumeStar (y = [x1 x2])
        function volume = concatenation(obj, I)
            
            S1 = obj.toStar;
            S2 = I.toStar;
            S = S1.concatenate(S2);

            volume = S.toVolumeStar(I.height, I.width, I.depth, I.numChannel + obj.numChannel);

        end

        % compute the hadamardProduct of two VolumeStars
        function volume = HadamardProduct(obj, I)
            % convert to Star
            S1 = obj.toStar;
            S2 = I.toStar;

            % Compute product
            S = S1.HadamardProduct(S2);

            % create new VolumeStar
            volume = S.toVolumeStar(I.height, I.width, I.depth, I.numChannel);
        end

        % transform to Star
        function S = toStar(obj)
            % get info about VolumeStar
            nc = obj.numChannel;
            h = obj.height;
            w = obj.width;
            dp = obj.depth;
            np = obj.numPred;
            
            N = h*w*dp*nc; % total number of pixels in the input volume        
            V1(:, np+1) = cast(zeros(N, 1), 'like', obj.V);
            for j=1:np+1
                V1(:,j) = reshape(obj.V(:, :, :, :, j), N, 1);
            end       
            if ~isempty(obj.vol_lb) && ~isempty(obj.vol_ub)
                state_lb = reshape(obj.vol_lb, N, 1);
                state_ub = reshape(obj.vol_ub, N, 1);
                S = Star(V1, obj.C, obj.d, obj.pred_lb, obj.pred_ub, state_lb, state_ub);
            else
                S = Star(V1, obj.C, obj.d, obj.pred_lb, obj.pred_ub);
            end
            
        end
        
        % checking if a VolumeStar is an empty set
        function bool = isEmptySet(obj)
            
            try
                S = obj.toStar;
                bool = S.isEmptySet;
            catch
                [lb, ub] = obj.getRanges;
                if isempty(lb) && isempty(ub)
                    bool = 1;
                else
                    bool = 0;
                end
            end
        end
        
        % contain, check if VolumeStar contains a volume image
        function bool = contains(obj, volume)
            % @volume: input volume image
            % @bool: = 1 if the volume is contained in the set
            %        = 0 if the volume is NOT contained in the set
            
            n = size(volume);
            if length(n) == 3 % one channel volume
                if n(1) ~= obj.height || n(2) ~= obj.width || n(3) ~= obj.depth || obj.numChannel ~= 1
                    error('Inconsistent dimenion between input volume and the VolumeStar');
                end
                y = reshape(volume, [n(1)*n(2)*n(3) 1]);
            elseif length(n) == 4
                if n(1) ~= obj.height || n(2) ~= obj.width || n(3) ~= obj.depth || n(4) ~= obj.numChannel
                    error('Inconsistent dimenion between input volume and the VolumeStar');
                end
                y = reshape(volume, [n(1)*n(2)*n(3)*n(4) 1]);
            else
                error('Invalid input volume');
            end
            
            % Convert to Star
            S = obj.toStar;
            % check for containment
            bool = S.contains(y);
            
        end
        
        % projection of volume on specific 2d plane
        function proj_star = project2D(obj, point1, point2)
            % @point1: index of first dimension
            % @point2: index of second dimension
            % @proj_star: is a projected star set
            
            if (length(point1) ~= 4) || (length(point2) ~= 4)
                error('Invalid input points. Points must by 4-D arrays');
            end
            
            if (point1(1) < 1) || (point1(1) > obj.height) ||...
                    (point1(2) < 1) || (point1(2) > obj.width) ||...
                    (point1(3) < 1) || (point1(3) > obj.depth) ||...
                    (point1(4) < 1) || (point1(4) > obj.numChannel)
                error('The first input point is invalid');
            end
            
            if (point2(1) < 1) || (point2(1) > obj.height) ||...
                    (point2(2) < 1) || (point2(2) > obj.width) ||...
                    (point2(3) < 1) || (point2(3) > obj.depth) ||...
                    (point2(4) < 1) || (point2(4) > obj.numChannel)
                error('The second input point is invalid');
            end
            
            n = obj.numPred + 1; 
            
            new_V = zeros(2, n);
            for i=1:n
                new_V(1, i) = obj.V(point1(1), point1(2), point1(3), point1(4), i);
                new_V(2, i) = obj.V(point2(1), point2(2), point2(3), point2(4), i);              
            end     
            
            proj_star = Star(new_V, obj.C, obj.d, obj.pred_lb, obj.pred_ub);
            
        end
        
        % TODO: add a 3D projection to ImageStar
    end


    methods % get methods
        
        % get ranges of a state at specific position
        function [xmin, xmax] = getRange(varargin)
            % [xmin, xmax] = getRange(vert_ind, horiz_ind, depth_ind, channel_ind)
            % INPUTS
            %   @vert_ind:  vectical index
            %   @horiz_ind: horizontal index
            %   @depth_ind: depth index
            %   @chan_ind : channel index
            % -----------------------------
            % OUTPUTS:
            %   @xmin: min of x(vert_ind, horiz_ind, depth_ind, channel_ind)
            %   @xmax: max of x(vert_ind, horiz_ind, depth_ind, channel_ind)
            
            switch nargin
                case 5
                    obj = varargin{1};
                    vert_ind = varargin{2};
                    horiz_ind = varargin{3};
                    depth_ind = varargin{4};
                    chan_ind = varargin{5};
                    lp_solver = 'linprog';
                case 6
                    obj = varargin{1};
                    vert_ind = varargin{2};
                    horiz_ind = varargin{3};
                    depth_ind = varargin{4};
                    chan_ind = varargin{5};
                    lp_solver = varargin{6};
                otherwise
                    error('Invalid number of input arguments, should be 3 or 4');
            end
            
            % Check input validity
            if isempty(obj.C) || isempty(obj.d)
                error('The VolumeStar is empty');
            end
            
            if vert_ind < 1 || vert_ind > obj.height
                error('Invalid veritical index');
            end
            
            if horiz_ind < 1 || horiz_ind > obj.width
                error('Invalid horizonal index');
            end

            if depth_ind < 1 || depth_ind > obj.depth
                error('Invalid depth index');
            end
            
            if chan_ind < 1 || chan_ind > obj.numChannel
                error('Invalid channel index');
            end

            % min
            f = obj.V(vert_ind, horiz_ind, depth_ind, chan_ind, 2:obj.numPred + 1);
            [fval, exitflag] = lpsolver(f, obj.C, obj.d, [], [], obj.pred_lb, obj.pred_ub, lp_solver);
            if ismember(exitflag, ["l1","g5"])
               xmin = fval + obj.V(vert_ind, horiz_ind, depth_ind, chan_ind, 1);
            else
                error("Cannot find an optimal solution, exitflag = " + string(exitflag));
            end

            % max
            [fval, exitflag] = lpsolver(-f, obj.C, obj.d, [], [], obj.pred_lb, obj.pred_ub, lp_solver);
            if ismember(exitflag, ["l1","g5"])
                xmax = -fval + obj.V(vert_ind, horiz_ind, depth_ind, chan_ind, 1);
            else
                error("Cannot find an optimal solution, exitflag = " + string(exitflag));
            end
            
            % update VolumeStar bounds
            obj.vol_lb(vert_ind, horiz_ind, depth_ind, chan_ind) = xmin;
            obj.vol_ub(vert_ind, horiz_ind, depth_ind, chan_ind) = xmax;
                   
        end
        
        % estimate range quickly using only predicate bound information
        function [xmin, xmax] = estimateRange(obj, h, w, dp, c)
            % @h: height index
            % @w: width index
            % @dp: depth index
            % @c: channel index
            % @xmin: min of x[h, w, dp, c]
            % @xmax: max of x[h, w, dp, c]
            
            if isempty(obj.C) || isempty(obj.d)
                error('The VolumeStar is empty');
            end
            
            if h < 1 || h > obj.height
                error('Invalid veritical index');
            end
            
            if w < 1 || w > obj.width
                error('Invalid horizonal index');
            end

            if dp < 1 || dp > obj.depth
                error('Invalid depth index');
            end
            
            if c < 1 || c > obj.numChannel
                error('Invalid channel index');
            end
            
            f = obj.V(h, w, dp, c, 1:obj.numPred + 1);
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
        function [volume_lb, volume_ub] = estimateRanges(varargin)
            % [volume_lb, volume_ub] = obj.estimateRanges;
            % @volume_lb: lower bound volume
            % @volume_ub: upper bound volume
            
            switch nargin
                case 1
                    obj = varargin{1};
                    dis_opt = [];
                case 2
                    obj = varargin{1};
                    dis_opt = varargin{2};
                otherwise
                    error('\nInvalid number of input arguments, should be 0 or 1');
            end
            
            if isempty(obj.C) || isempty(obj.d)
                error('The VolumeStar is empty');
            end
            
            if isempty(obj.vol_lb) || isempty(obj.vol_ub)
                % initialize vars
                volume_lb = zeros(obj.height, obj.width, obj.depth, obj.numChannel);
                volume_ub = zeros(obj.height, obj.width, obj.depth, obj.numChannel);
                reverseStr = '';
                % get total number of pixels
                N = obj.height*obj.width*obj.depth*obj.numChannel;
                if strcmp(dis_opt, 'display')
                    fprintf('\nEstimating Range: ');
                end
                % Estimate range pixel by pixel
                for h = 1:obj.height
                    for w = 1:obj.width
                        for dp = 1:obj.depth
                            for c=1:obj.numChannel
                                [volume_lb(h, w, dp, c), volume_ub(h, w, dp, c)] = obj.estimateRange(h, w, dp, c);
                                if strcmp(dis_opt, 'display')
                                    msg = sprintf('%d/%d', h*w*dp*c, N);   
                                    fprintf([reverseStr, msg]);
                                    reverseStr = repmat(sprintf('\b'), 1, length(msg));
                                end
                            end
                        end
                    end
                end
                % update bounds
                obj.vol_lb = volume_lb;
                obj.vol_ub = volume_ub;
            else
                % get precomputed bunds
                volume_lb = obj.vol_lb;
                volume_ub = obj.vol_ub;
            end
         
        end
        
        % get lowew bound and upper bound volumes of an VolumeStar
        function [volume_lb, volume_ub] = getRanges(varargin)
            % [volume_lb, volume_ub] = obj.getRanges;
            % @volume_lb: lower bound volume
            % @volume_ub: upper bound volume
            
            switch nargin
                case 1
                    obj = varargin{1};
                    lp_solver = 'linprog';
                case 2
                    obj = varargin{1};
                    lp_solver = varargin{2};
                otherwise
                    error('Invalid number of input arguments');
            end
            
            volume_lb = zeros(obj.height, obj.width, obj.depth, obj.numChannel);
            volume_ub = zeros(obj.height, obj.width, obj.depth, obj.numChannel);

            % Estimate range pixel by pixel
                for h = 1:obj.height
                    for w = 1:obj.width
                        for dp = 1:obj.depth
                            for c=1:obj.numChannel
                                [volume_lb(h, w, dp, c), volume_ub(h, w, dp, c)] = obj.getRange(h, w, dp, c, lp_solver);
                            end
                        end
                    end
                end
            
            obj.vol_lb = volume_lb;
            obj.vol_ub = volume_ub;
                   
        end
        
        % estimate the number of attacked pixels
        function n_att = getNumAttackedPixels(obj)
            % @n_att: number of attacked pixels in an VolumeStar
            
            V1 = zeros(obj.height, obj.width, obj.depth, obj.numChannel);
            V3 = V1;
            for i=2:obj.numPred + 1
                V2 = (obj.V(:,:,:,i) ~= V1);
                V3 = V3 + V2; 
            end          
            [V4, ~] = max(V3, [], 4);
            n_att = sum(V4, 'all');
        end
        
    end


    methods % Pooling related operation
        
        % update local ranges for Max Pooling operation
        function updateRanges(varargin)
            % @points: local points = [x1 y1 z1 c1; x2 y2 z2 c2; ...]
            
            switch nargin
                case 2
                    obj = varargin{1};
                    points = varargin{2};
                    lp_solver = 'linprog';
                case 3
                    obj = varargin{1};
                    points = varargin{2};
                    lp_solver = varargin{3};
                otherwise
                    error('Invalid number of input arguments, should be 2 or 3');
            end

            n = size(points, 1);
            for i=1:n
                obj.getRange(points(i, 1), points(i,2), points(i, 3), points(i,4), lp_solver);
            end        
            
        end

        % get local bound for Max Pooling operation
        function [lb, ub] = get_localBound(varargin)
            % @startpoint: startpoint of the local volume
            %               startpoint = [x1 y1];
            % @PoolSize: = [height width depth]  of max pooling layer
            % @channel_id: the index of the channel
            % @lb, ub: the lower bound and upper bound of all points in the local region
            
            switch nargin
                case 4
                    obj = varargin{1};
                    startpoint = varargin{2};
                    PoolSize = varargin{3};
                    channel_id = varargin{4};
                    lp_solver = 'linprog';
                case 5
                    obj = varargin{1};
                    startpoint = varargin{2};
                    PoolSize = varargin{3};
                    channel_id = varargin{4};
                    lp_solver = varargin{5};
                otherwise
                    error('Invalid number of input arguments, should be 4 or 5');
            end
            
            points = obj.get_localPoints(startpoint, PoolSize);
            n = length(points);
            % get lower bound and upper bound volume
            if isempty(obj.vol_lb) || isempty(obj.vol_ub)
                [volume_lb, volume_ub] = obj.getRanges(lp_solver);
            else
                volume_lb = obj.vol_lb;
                volume_ub = obj.vol_ub;
            end
            
            lb = volume_lb(points(1,1), points(1,2), points(1,3), channel_id);
            ub = volume_ub(points(1,1), points(1,2), points(1,3), channel_id);
            
            for i=2:n
                if volume_lb(points(i,1), points(i,2), points(i,3), channel_id) < lb
                    lb = volume_lb(points(i,1), points(i,2), points(i,3), channel_id);
                end
                if volume_ub(points(i,1), points(i,2), points(i,3), channel_id) > ub
                    ub = volume_ub(points(i,1), points(i,2), points(i,3), channel_id);
                end
            end
            
        end
        
        % get all local points index for Max Pooling operation
        function points = get_localPoints(obj, startpoint, PoolSize) % TODO: check permuted array is correctly reshaped (dimension order)
            % @startpoint: startpoint of the local volume (startpoint = [x1 y1 z1]);
            % @PoolSize: = [height width] the height and width of max pooling layer
            % @points: all indexes of all points for a single max
            % pooling operation (including the startpoint)
            
            x0 = startpoint(1); % vertical index of the startpoint
            y0 = startpoint(2); % horizontal index of the startpoint
            z0 = startpoint(3); % depth index of the startpoint
            h  = PoolSize(1);   % height of the MaxPooling layer
            w  = PoolSize(2);   % width of the MaxPooling layer
            dp = PoolSize(3);   % depth of the MaxPooling layer
            
            if x0 < 1 || y0 < 1 || z0 < 1 || x0 + h - 1 > obj.height || ...
                    y0 + w - 1 > obj.width || z0 + dp - 1 > obj.depth
                error('Invalid startpoint or PoolSize');
            end
            points = zeros(h*w*dp);
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

                    for k=1:dp
                        if j==1
                        z1 = z0;
                        else
                            z1 = z1 + 1;
                        end

                        points(i,j,k) = [x1 y1 z1];
                    end
                end
            end
            points = reshape(points, [], h*w*dp, 3);
        end
            
        % get local max index( find the maximum point of a local volume)
        function max_id = get_localMax_index(varargin) %% TODO: begin here again
            % @startpoint: startpoint of the local volume (startpoint = [x1 y1 z1]);
            % @PoolSize: = [height width depth] of max pooling layer
            % @channel_id: the channel index
            % @max_id: = []: we don't know which one has maximum value,
            % i.e., the maximum values may be the intersection between of
            % several pixel values => [xi yi xi]: the point that has maximum value
            % (used in over-approximate reachability analysis of maxpooling operation)
            
            switch nargin
                case 4
                    obj = varargin{1};
                    startpoint = varargin{2};
                    PoolSize = varargin{3};
                    channel_id = varargin{4};
                    lp_solver = 'linprog';
                case 5
                    obj = varargin{1};
                    startpoint = varargin{2};
                    PoolSize = varargin{3};
                    channel_id = varargin{4};
                    lp_solver = varargin{5};
                otherwise
                    error('Invalid number of input arguments, should be 4 or 5');
            end
            
            points = obj.get_localPoints(startpoint, PoolSize);          
            % get lower bound and upper bound volume
            if isempty(obj.vol_lb) || isempty(obj.vol_ub)
                obj.estimateRanges;
            end
            
            n = prod(PoolSize); % number of pixels (flattened to vector)
            
            lb = zeros(1, n);
            ub = zeros(1, n);
            % get ubber and lower bound vectors
            for i=1:n
                point_i = points(i, :);
                lb(:, i) = obj.vol_lb(point_i(1), point_i(2), point_i(3), channel_id);
                ub(:, i) = obj.vol_ub(point_i(1), point_i(2), point_i(3), channel_id);   
            end
            % get max values
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
                obj.updateRanges(new_points, lp_solver);
                                
                lb = zeros(1, m);
                ub = zeros(1, m);

                for i=1:m
                    point_i = points(candidates(i), :);
                    lb(:, i) = obj.vol_lb(point_i(1), point_i(2), point_i(3), channel_id);
                    ub(:, i) = obj.vol_ub(point_i(1), point_i(2), point_i(3), channel_id);   
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
                        if obj.is_p1_larger_p2([p1(1) p1(2) p1(3) channel_id], [max_id(1) max_id(2) max_id(3) channel_id], lp_solver)
                            max_id1 = [max_id1; p1];
                        end
                    end                
                    max_id = max_id1; 
                    
                end
                
            end
            
            n = size(max_id, 1);
            channels = channel_id*ones(n,1);
            max_id = [max_id channels];

        
        end
        
        % get local max index (find the maximum point of a local volume) 
        function max_id = get_localMax_index2(obj, startpoint, PoolSize, channel_id)
            % @startpoint: startpoint of the local volume (startpoint = [x1 y1 z1])
            % @PoolSize: = [height width depth] of max pooling layer
            % @channel_id: the channel index
            % @max_id: = []: we don't know which one has maximum value,
            % i.e., the maximum values may be the intersection between of
            % several pixel valutes => [xi yi zi]: the point that has maximum value
            % used in over-approximate reachability analysis of maxpooling operation
            
            points = obj.get_localPoints(startpoint, PoolSize);          
            % get lower bound and upper bound volume
            if isempty(obj.vol_lb) || isempty(obj.vol_ub)
                obj.estimateRanges;
            end
            
            n = prod(PoolSize);
            
            lb = zeros(1, n);
            ub = zeros(1,n);
            
            for i=1:n
                point_i = points(i, :);
                lb(:, i) = obj.vol_lb(point_i(1), point_i(2), point_i(3), channel_id);
                ub(:, i) = obj.vol_ub(point_i(1), point_i(2), point_i(3), channel_id);   
            end
            
            [max_lb_val, ~] = max(lb, [], 2);
            
            max_id = find(ub - max_lb_val > 0);
            n = size(max_id, 1);
            channels = channel_id*ones(n,1);
            max_id = [max_id channels];
      
        end
        
        % add maxidx (used for unmaxpooling reachability)
        function addMaxIdx(obj, name, maxIdx)
            % @name: name of the max pooling layer
            % @maxIdx: max indexes
            
            A.Name = name;
            A.MaxIdx = maxIdx;
            if isempty(obj.MaxIdxs{1})
                obj.MaxIdxs{1} = A;
            else
                obj.MaxIdxs = [obj.MaxIdxs A];
            end
        end
        
        % update max index (used for unmaxpooling reachability)
        function updateMaxIdx(obj, name, maxIdx, pos)
            % @name: name of the max pooling layer
            % @maxIdx: max indexes
            % @pos: the position of the local pixel of the max map
            % corresponding to the maxIdx
            
            n = length(obj.MaxIdxs);
            ct = 0;
            for i=1:n
                if strcmp(obj.MaxIdxs{i}.Name, name)
                    obj.MaxIdxs{i}.MaxIdx{pos(1), pos(2), pos(3), pos(4)} = maxIdx;
                    break;
                else
                    ct = ct + 1;
                end
            end
            
            if ct == n
                error('Unknown name of the maxpooling layer');
            end
            
        end
        
        % add input size (used for unmaxpooling reachability)
        function addInputSize(obj, name, inputSize)
            % @name: name of the max pooling layer
            % @inputSize: input size of the original volume
                       
            A.Name = name;
            A.InputSize = inputSize;
            if isempty(obj.InputSizes{1})
                obj.InputSizes{1} = A;
            else
                obj.InputSizes = [obj.InputSizes A];
            end
        end
        
        % compare between two specific points in the volume (is p1 > p2 isfeasible?)
        function b = is_p1_larger_p2(varargin)
            % @p1: the first point  = [h1, w1, dp1, c1]
            % @p2: the second point = [h2, w2, dp2, c2]
            % h: height, w: width, dp: depth, c: channel index
            %
            % @b = 1 -> p1 > p2 is feasible
            %    = 0 -> p1 > p2 is not feasible
            % % useful for max pooling operation
            
            switch nargin
                case 3
                    obj = varargin{1};
                    p1 = varargin{2};
                    p2 = varargin{3};
                    lp_solver = 'linprog';
                case 4
                    obj = varargin{1};
                    p1 = varargin{2};
                    p2 = varargin{3};
                    lp_solver = varargin{4};
                otherwise
                    error('Invalid number of input arguments, should be 2 or 3');
            end
            
            C1 = zeros(1, obj.numPred);
            for i=2:obj.numPred+1
                C1(i-1) = obj.V(p2(1), p2(2), p2(3), p2(4), i) - obj.V(p1(1), p1(2), p1(3), p1(4), i);
            end
            
            d1 = obj.V(p1(1), p1(2), p1(3), p1(4), 1) - obj.V(p2(1), p2(2), p2(3), p2(4)); 
            
            new_C = [obj.C; C1];
            new_d = [obj.d; d1];
           
            f = zeros(1, obj.numPred);
            [~, exitflag] = lpsolver(f, new_C, new_d, [], [], obj.pred_lb, obj.pred_ub, lp_solver);
            if ismember(exitflag, ["l1", "g2", "g5"]) % feasible solution exist
                b = 1;
            elseif ismember(exitflag, ["l-2", "l-5", "g3", "g4", "g110"])
                b = 0;
            else
                error("ERROR, exitflag = " + string(exitflag));
            end
            
        end
         
    end
    
    
    methods(Static) % helper function
        
        % check if a pixel value is the maximum value compared with others
        function [new_C, new_d] = isMax(varargin)
            % @maxMap: the current maxMap VolumeStar
            % @ori_volume: the original VolumeStar to compute the maxMap 
            % @center: is the center pixel position we want to check
            %          center = [x1 y1 z1 c1]
            % @others: is the other pixel position we want to compare with
            % the center one
            %          others = [x2 y2 z2 c2; x3 y3 z3 c3]
            %
            % this is core step for exactly performing maxpooling operation on an VolumeStar set
            
            switch nargin
                case 4 
                    maxMap = varargin{1};
                    ori_volume = varargin{2};
                    center = varargin{3};
                    others = varargin{4};
                    lp_solver = 'linprog';
                case 5
                    maxMap = varargin{1};
                    ori_volume = varargin{2};
                    center = varargin{3};
                    others = varargin{4};
                    lp_solver = varargin{5};
                otherwise 
                    error('Invalid number of input arguments, should be 4 or 5');
            end
                       
            if maxMap.numPred ~= ori_volume.numPred
                error('Inconsistency between number of predicates in the current maxMap and the original volume');
            end
            
            n = size(others, 1);
           
            % the center may be the max point with some extra constraints on the predicate variables
            new_C = zeros(n, maxMap.numPred);
            new_d = zeros(n, 1);
            
            for i=1:n                
                % add new constraint
                new_d(i) = ori_volume.V(center(1),center(2), center(3), center(4), 1) - ori_volume.V(others(i,1), others(i,2), others(i,3), others(i,4), 1);
                for j=1:maxMap.numPred                    
                    new_C(i,j) = -ori_volume.V(center(1),center(2), center(3), center(4), j+1) + ori_volume.V(others(i,1), others(i,2), others(i,3), others(i,4), j+1);
                end           
            end
            
            C1 = [maxMap.C; new_C];
            d1 = [maxMap.d; new_d];

            % remove redundant constraints
            E = [C1 d1];
            E = unique(E, 'rows');

            C1 = E(:, 1:ori_volume.numPred);
            d1 = E(:, ori_volume.numPred + 1);

            f = zeros(1, ori_volume.numPred);
            [~, exitflag] = lpsolver(f, C1,  d1, [], [], ori_volume.pred_lb, ori_volume.pred_ub, lp_solver);
            if ismember(exitflag, ["l1", "g5"]) % optimal solution exist
                new_C = C1;
                new_d = d1;
            else
                new_C = [];
                new_d = [];
            end

        end
                        
        % reshape an VolumeStar
        function new_VS = reshape(inputSet, new_shape)
            % @inutSet: input VolumeStar
            % @new_shape: new shape
            % @new_VS: new VolumeStar

            n = size(new_shape);

            if n(2) ~= 4 || n(1) ~= 1
                error('Target shape should be 1x4 array.');
            end

            if prod(new_shape) ~= inputSet.height * inputSet.width * inputSet.depth * inputSet.numChannel
                error('New shape is inconsistent with the current shape');
            end
            
            new_V = reshape(inputSet.V, [new_shape inputSet.numPred + 1]);           
            new_VS = VolumeStar(new_V, inputSet.C, inputSet.d, inputSet.pred_lb, inputSet.pred_ub, inputSet.vol_lb, inputSet.vol_ub);
     
        end
        
        % add new constraint to predicate variables of an VolumeStar
        function [new_C, new_d] = addConstraint(inputSet, p1, p2)
            % @inputSet: input VolumeStar
            % @p1: first point position
            % @p2: second point position
            % @new_C: a new predicate constraint matrix C 
            % @new_d: new predicate constraint vector
            %
            % used for finding counter examples
            % we will add new constraint: p2 >= p1
            
            if ~isa(inputSet, 'VolumeStar')
                error('Input set is not an VolumeStar');
            end
            
            new_d = inputSet.V(p2(1), p2(2), p2(3), p2(4), 1) - inputSet.V(p1(1), p1(2), p1(3), p1(4), 1);
            new_C = inputSet.V(p1(1), p1(2), p1(3), p1(4), 2:inputSet.numPred + 1) - inputSet.V(p2(1), p2(2), p2(3), p2(4), 2:inputSet.numPred + 1);
            
            new_C = reshape(new_C, [1 inputSet.numPred]);
            new_C = [inputSet.C; new_C];
            new_d = [inputSet.d; new_d];
            
        end

    end
    
    
end


