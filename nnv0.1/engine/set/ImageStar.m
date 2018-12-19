classdef ImageStar
    % Class for representing set of images using Star set
    % An image can be attacked by bounded noise. An attacked image can
    % be represented using an ImageStar Set
    % Dung Tran: 12/17/2018
    
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
        
        % 2D representation of an ImageStar
        IM = []; % center image (high-dimensional array)
        LB = []; % lower bound of attack (high-dimensional array)
        UB = []; % upper bound of attack (high-dimensional array)
        
        % 1D representation of an ImageStar
        % flattening image to a normal star set: S = c + V*a (see Star class)
        Stars = []; % an array of stars, number of stars = numChannel
        

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
                case 3
                    IM = varargin{1};
                    LB = varargin{2};
                    UB = varargin{3};
                    n = size(IM); % n(1) and n(2) are height and width of image
                                  % n(3) is number of channels
                    l = size(LB);
                    u = size(UB);

                    if n(1) ~= l(1) || n(1) ~= u(1) || n(2) ~= l(2) || n(2) ~= u(2) 
                        error('Inconsistency between center image and attack bound matrices');
                    end

                    if length(n) ~= length(l) || length(n) ~= length(u)
                        error('Inconsistency between center image and attack bound matrices');
                    end

                    if length(n) == 2 && length(l) == 2 && length(u) == 2

                        obj.numChannel = 1;
                        obj.IM = IM;
                        obj.LB = LB;
                        obj.UB = UB;
                        obj.height = n(1);
                        obj.width = n(2);

                    elseif length(n) == 3 && length(l) == 3 && length(u) == 3

                        if n(3) == l(3) && n(3) == u(3)
                            obj.numChannel = n(3);
                            obj.IM = IM; 
                            obj.LB = LB;
                            obj.UB = UB;
                            obj.height = n(1);
                            obj.width = n(2);
                        else
                            error('Inconsistent number of channels between the center image and the bound matrices');
                        end

                    else
                        error('Inconsistent number of channels between the center image and the bound matrices');
                    end

                    % flattening 2D ImageStar to get 1D ImageStar (a normal star set)
                    S = [];
                    for i=1:obj.numChannel
                        c = reshape(obj.IM(:,:,i)', [obj.height * obj.width,1]);
                        lb = reshape(obj.LB(:,:,i)', [obj.height * obj.width,1]);
                        ub = reshape(obj.UB(:,:,i)', [obj.height * obj.width,1]);
                        lb = lb + c;
                        ub = ub + c;
                        B = Box(lb, ub);
                        S = [S B.toStar];

                    end
                    obj.Stars = S;
                    
                case 2
                    
                    S = varargin{1}; % 1D representation of an ImageStar
                    imageSize = varargin{2}; % size of the image
                    
                    if length(imageSize) ~= 2
                        error('Invalid image size of an ImageStar set');
                    end
                    
                    obj.height = imageSize(1);
                    obj.width = imageSize(2);
                    obj.numChannel = length(S);
                    obj.Stars = S; 
                    
                    for i=1:obj.numChannel
                        if obj.Stars(i).dim ~= obj.height * obj.width
                            error('Inconsistency between dimension of Stars and imageSize');
                        end
                    end
                    
                                 
                otherwise
                    
                    % create an empty ImageStar
                    obj.numChannel = 0; 
                    obj.height = 0;
                    obj.width = 0;
                    obj.IM = [];
                    obj.LB = [];
                    obj.UB = [];
                    obj.Stars = []; 
                    
            end
                        
        end
        
        
        % zero-padding an ImageStar
        % zero-padding an ImageStar results another ImageStar
        function padded_image = zero_padding(obj, paddingSize)
            % @paddingSize: [t b l r] an 1-D array
            % @image: a new imageStar after padding           
            % @reference: https://www.mathworks.com/help/deeplearning/ug/layers-of-a-convolutional-neural-network.html
            
            % author: Dung Tran
            % date: 12/17/2018
            
            n = size(paddingSize);
            if n(1) ~= 1
                error('Padding Size should be a one row matrix');
            end
            if n(2) ~= 4
                error('Padding Size should have 4 column');
            end

            t = paddingSize(1); % top padding
            b = paddingSize(2); % botton padding
            l = paddingSize(3); % left padding
            r = paddingSize(4); % right padding

            if t < 0 || b < 0 || l < 0 || r < 0
                error('Invalid padding size');
            end
            
            % zero-padding for 2D representation of an image star
            if ~isempty(obj.IM) && ~isempty(obj.LB) && ~isempty(obj.UB) 
                h = obj.height + t + b; % height of new image
                w = obj.width + l + r; % width of new image
                new_IM = zeros(h, w, obj.numChannel); % preallocate new image
                new_LB = zeros(h, w, obj.numChannel); % preallocate new lower bound matrix
                new_UB = zeros(h, w, obj.numChannel); % preallocate new upper bound matrix
                for i=1:obj.numChannel
                    new_IM(t+1:t+obj.height, l+1:l+obj.width, i) = obj.IM(:,:,i);
                    new_LB(t+1:t+obj.height, l+1:l+obj.width, i) = obj.LB(:,:,i);
                    new_UB(t+1:t+obj.height, l+1:l+obj.width, i) = obj.UB(:,:,i);
                end

                padded_image = ImageStar(new_IM, new_LB, new_UB);

            % zero-padding for 1D representation of an image star  
            elseif ~isempty(obj.Stars)
                
                n = length(obj.Stars);
                 
                new_Stars(n) = obj.Stars(1); % preallocating stars array (faster performance without preallocation)
                for i=1:n
                   
                    V = obj.Stars(i).V;
                    m = size(V, 2);
                    
                    m1 = (t + b + obj.height) * (l + r + obj.width);
                    V1 = zeros(m1, m);
                    
                    for j=1:m
                        C1 = zeros(t + b + obj.height, l + r + obj.width); 
                        C = V(:, j);
                        C = reshape(C, [obj.height, obj.width]);
                        C1(t+1:t+obj.height, l+1:l + obj.width) = C';
                        V1(:, j) = reshape(C1', [m1, 1]);                       
                    end
                    new_Stars(i) = Star(V1, obj.Stars(i).C, obj.Stars(i).d);
                    
                end
                
                padded_image = ImageStar(new_Stars, [t + b + obj.height, l + r + obj.width]);
                
                
            else
                error('The image star is empty');
            end
            
            
        end
        
        
        % Filtering an ImageStar at specific region by a filter W
        % Reference: https://www.mathworks.com/help/deeplearning/ug/layers-of-a-convolutional-neural-network.html
        % =================================================================%
        % ***Importance: A filtered Image Star is an ImageStar that does
        % not has a full 2D representation. LB = [] and UB = []
        % The reason is that: 2D representation can not precisely represent
        % a filtered ImageStar. To precisely reprent a filtered ImageStar,
        % we use 1D representation. 
        % we can obtain a 2D presentation from 1D representation. However,
        % in this case, 2D representation is an over-approximation (box) of
        % the exact 1D representation. 
        
        % We always use 1D representation to compute the exact reachable set 
        % of a convolutional 2D layer. 
        % =================================================================%
        
        
        
        
               
    end
end

