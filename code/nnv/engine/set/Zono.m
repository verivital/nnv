classdef Zono
    %Zonotope class 
    %   Z = (c , <v1, v2, ..., vn>) = c + a1 * v1 + ... + an * vn, 
    %   where -1 <= ai <= 1
    %   c is center, vi is generator
    % Dung Tran: 10/23/2018
    
    properties
        c = []; % center vector
        V = []; % generator matrix V = [v1 v2 ... vn]
        dim = 0; % dimension of the zonotope
    end
    
    methods % constructor and main methods
        
        %constructor
        function obj = Zono(varargin)
            % @c: center vector
            % @V: generator matrix
            
            switch nargin
                
                case 2
                    c = varargin{1};
                    V = varargin{2};
                    [nC, mC] = size(c);
                    [nV, ~] = size(V);

                    if mC ~= 1
                        error('center vector should have one column');
                    end
                    if nC ~= nV
                        error('Inconsistent dimension between center vector and generator matrix');
                    end

                    obj.c = c;
                    obj.V = V;
                    obj.dim = nV;
                case 0 
                otherwise
                    error('Invalid number of inputs, 0 or 2');
            end
            
        end
        
        % affine mapping of a zonotope: Wz + b     
        function Z = affineMap(obj, W, b)
            % @W: affine mapping matrix
            % @b: mapping vector
            
            [nW, mW] = size(W);
            [nb, mb] = size(b);
            
            if mW ~= obj.dim
                error('Inconsistent dimension between mapping matrix with the zonotope dimension');
            end
            
            if ~isempty(b)
                if mb ~= 1
                    error('Mapping vector should have one column');
                end
                if nb ~= nW
                    error('Inconsistent dimension between mapping vector and mapping matrix');
                end
                new_c = W * obj.c + b;
                new_V = W * obj.V;
            else
                new_c = W * obj.c;
                new_V = W * obj.V;
            end

            Z = Zono(new_c, new_V);
            
        end
        
        % Minkowski Sum of two zonotopes
        function Z = MinkowskiSum(obj, X)
            % @X: input zonotope
            % @Z: output zonotope 
            
            if ~isa(X, 'Zono')
                error('Input set is not a zonotope');
            else
                if X.dim ~= obj.dim
                    error('Inconsistent dimension of input zonotope and this zonotope');
                end
            end
            
            Z = Zono(obj.c + X.c, horzcat(obj.V, X.V));
            
        end
        
        % convex hull with another zonotope
        function Z = convexHull(obj, X)
            % @X: input zonotope
            % @Z: output zonotope
            
            % author: Dung Tran
            % date: 10/25/2018
            % main reference: 1) Reachability of Uncertain Linear Systems Using
            % Zonotopes, Antoin Girard, HSCC2005.
            
            % ===================================================================== #
            % Convex hull of two zonotopes is generally NOT a zonotope.
            % This method returns an over-approximation (which is a zonotope)
            % of a convex hull of two zonotopes.
            % This method is a generalization of the method proposed in reference 1.
            % In reference 1: the author deals with: CONVEX_HULL(Z, L*Z)
            % Here we deals with a more general case: CONVEX_HULL(Z1, Z2)
            % We will see that if Z2 = L * Z1, the result is reduced to the result
            % obtained by the reference 1.
            % We define CH as a convex hull operator and U is union operator.
            % Z1 = (c1, <g1, g2, ..., gp>), Z2 = (c2, <h1, ...., hq>)
            % ===================================================================== #
            % CH(Z1 U Z2) := {a * x1 + (1 - a) * x2}| x1 \in Z1, x2 \in Z2, 0 <= a <= 1}
            % Let a = (e + 1)/2, -1<= e <=1, we have:
            %     CH(Z1 U Z2) := {(x1 + x2)/2 + e * (x1 - x2)/2}
            %                  = (Z1 + Z2)/2 + e*(Z1 + (-Z2))/2
            %                  where, '+' denote minkowski sum of two zonotopes
            % From minkowski sum method, one can see that:
            %      CH(Z1 U Z2) = 0.5*(c1 + c2 <2g1, ..., 2gp, 2h1, ...2hq, c1 - c2>)
            % So, the zonotope that over-approximate the convex hull of two zonotop
            % has (p + q) + 1 generators.
            %                                               
            % So, the zonotope that over-approximate the convex hull of two zonotop
            % has (p + q) + 1 generators.
            % Let consider the specific case Z2 = L * Z1.
            % In this case we have:
            %     (Z1 + L*Z1)/2 = 0.5 * (I+L) * (c1, <g1, g2, ..., gp>)
            %     (Z1 - L*Z1)/2 = 0.5 * (I-L) * (c1, <g1, ..., gp>)
            %     ea * (Z1 - L*Z1)/2 = 0.5*(I-L)*(0, <c1, g1, ..., gp>)
            %     CH(Z1 U L * Z1) = 0.5*((I + L)*c1, <(I+L)*g1, ..., (I+L)*gp,
            %                        (I-L)*c1, (I-L)*g1, ..., (I-L)*gp>)
            % where I is an identity matrix.
            % So the resulted zonotope has 2p + 1 generators.
            % ===================================================================== #

            if ~isa(X, 'Zono')
                error('Input set is not a zonotope');
            end
            
            if X.dim ~= obj.dim
                error('Inconsistent dimension between input set and this zonotope');
            end
            
            new_c = 0.5*(obj.c + X.c);
            new_V = [obj.V X.V 0.5*(obj.c -X.c)];
            Z = Zono(new_c, new_V);
            
        end
        
        % convexHull of a zonotope with its linear transformation
        function Z = convexHull_with_linearTransform(obj, L)
            % @L: linear transformation matrix
            % @Z: new zonotope
            
            % author: Dung Tran
            % date: 10/27/2018
            
            [nL, mL] = size(L);
            
            if nL ~= mL 
                error('Trasformation matrix should be a square matrix');
            end
            if nL ~= obj.dim
                error('Inconsistent dimension of tranformation matrix and this zonotope');
            end
            
            M1 = eye(nL) + L;
            M2 = eye(nL) - L;
            
            new_c = 0.5 * M1 * obj.c;
            new_V = 0.5*[M1 * obj.V M2 * obj.c M2 * obj.V];
            
            Z = Zono(new_c, new_V);
                        
        end
        
        % intersect with HalfSpace
        function S = intersectHalfSpace(obj, H, g)
            % @H: half space matrix
            % @g: half space vector
            % @S: intersection is a star
            
            % author: Dung Tran
            % date: 10/20/2019
            
            S = obj.toStar;
            S = S.intersectHalfSpace(H, g);            
        end

    end


    methods % conversion methods

        % Zonotope order reduction 
        function Z = orderReduction_box(obj, n_max)
            % @n_max: maximum allowable number of generators
            % @Z: new zonotope with order of n_max/dim
            
            % author: Dung Tran
            % date: 10/29/2018

            % Zonotope order reduction 
            % References 1) Reachability of Uncertain Linear Systems Using Zonotopes, Antoine Girard, HSCC 2008 (main reference)
            %            2) Methods for order reduction of Zonotopes, Matthias Althoff, CDC 2017
            %            3) A State bounding observer based on zonotopes, C. Combastel, ECC 2003
            % ====================================================================================== #
            %             ZONOTOPE ORDER REDUCTION USING BOX METHOD
            %
            % We define: R(Z) is the reduced-order zonotope that over-approximates
            % zonotope Z, i.e., Z \subset R(Z).
            %
            % IDEA: R(Z) is generated by less segments than Z.
            %
            % PRINCIPLE OF REDUCTION: "The edge of Z with lower length having priority
            % to be involved in the reduction". Read references 1, 3, 4, 5 for more information.
            %
            % STEPS:
            %    0) Z = (c, <g1, g2, ..., gp>) = c + Z0, Z0 = (0, <g1, g2, ..., gp>)
            %       we do order reduction for zonotope Z0, and then shift it to a new center point c.
            %
            %    1) Sort the generators g1, g2, ..., gp by their length to have:
            %       ||g1|| <= ||g2|| <= ||g3||.... <= ||gp||, where ||.|| is the 2-norm.
            %       This is the heuristic used in reference 3.
            %       we can use another heuristic used in reference 1 that is:
            %          ||g1||_1 - ||g1||_\inf <= .... <= ||gp||_1 - ||gp||_\inf.
            %       where ||.||_1 is 1-norm and ||.||_\inf is the infinity norm.
            %
            %    2) cur_order = p/n is the current order in which p is the number of generators,
            %       n is the dimension. r is the desired order of the reduced zonotope R(Z).
            %       r is usually selected as 1.
            %       Let d = cur_order - r be the number of orders needed to be reduced.
            %       We have p = (r + d)*n = (r-1)*n + (d+1)*n
            %
            %    3) The zonotope Z0 is splited into 2 zonotopes: Z0 = Z01 + Z02
            %       Z01 contains (d + 1)*n smallest generators
            %       Z02 contains (r - 1)*n lagest generators
            %
            %    4) Over-approximate Z01 by an interval hull IH(Z01), see references 3, 4, 5
            %       The IH(Z01) has n generators
            %
            %    5) The reduced zonotope is the Minkowski sum of IH(Z01) and Z02:
            %                    R(Z0) = IH(Z01) + Z02
            %       R(Z0) has (r-1)*n + n = r*n generators. So it has order of r.
            % ====================================================================================== #
            
            if n_max < obj.dim
                error('n_max should be >= %d', obj.dim);
            end      
            
            n = size(obj.V, 2); % number of generators
            
            if n <= n_max
                
                Z = Zono(obj.c, obj.V);
                
            else
                
                if n_max > obj.dim 
                    n1 = n - n_max + obj.dim; % number of generators need to be reduced
                   
                    % sort generators based on their lengths
                    length_gens = zeros(n, 1);
                    for i=1:n
                        length_gens(i) = norm(obj.V(:, i), 2);
                    end
                    [~, sorted_ind] = sort(length_gens); 

                    sorted_gens = zeros(obj.dim, n);
                    for i=1:n
                        sorted_gens(:, i) = obj.V(:, sorted_ind(i));
                    end

                    Z1 = Zono(zeros(obj.dim, 1), sorted_gens(:, 1:n1));
                    Z2 = Zono(obj.c, sorted_gens(:, n1+1:n));                          
                    Z = Z2.MinkowskiSum(Z1.getIntervalHull);
                   
                end
                
                if n_max == obj.dim
                    Z = obj.getIntervalHull();
                end
                
            end            
            
            
        end
        
        % convert to polyhedron
        function P = toPolyhedron(obj)
            n  = size(obj.V, 2);           
            C  = cast([eye(n); -eye(n)], 'like', obj.V);
            d  = cast(ones(2*n, 1), 'like', obj.V);
            Pa = Polyhedron(C, d);
            P  = Pa.affineMap(obj.V) + obj.c;
        end
        
        % convert to Star
        function S = toStar(obj)
            n = size(obj.V, 2);           
            lb = cast(-ones(n, 1), 'like', obj.V);
            ub = cast(ones(n, 1), 'like', obj.V);
            
            C = cast([eye(n);-eye(n)], 'like', obj.V);
            d = cast([ones(n,1); ones(n,1)], 'like', obj.V);
            S = Star([obj.c obj.V], C, d, lb, ub, obj);            
        end
        
        % convert to ImageZono
        function imageZono = toImageZono(obj, height, width, numChannels)
            % @height: height of the image
            % @width: width of the image
            % @numChannels: number of channels of the image
            % @imageZono: returned image
            
            % author: Dung Tran
            % date: 1/2/2020
            
            if height*width*numChannels ~= obj.dim
                error('Inconsistent dimension, please change the height, width and numChannels to be consistent with the dimension of the zonotope');
            end
             
            new_V = [obj.c obj.V]; 
            numPreds = size(obj.V, 2); 
            new_V = reshape(new_V, [height width numChannels numPreds + 1]);
            imageZono = ImageZono(new_V);
            
        end
        
        % convert to ImageStar
        function imageStar = toImageStar(obj, height, width, numChannels)
            % @height: height of the image
            % @width: width of the image
            % @numChannels: number of channels of the image
            % @imageZono: returned image
            
            % author: Dung Tran
            % date: 1/2/2020
            
            im1 = obj.toStar;
            imageStar = im1.toImageStar(height, width, numChannels);
            
        end
        
        % change variable precision
        function S = changeVarsPrecision(obj, precision)
            S = obj;
            if strcmp(precision, 'single')
                S.V = single(S.V);
                S.c = single(S.c);
            elseif strcmp(precision, 'double')
                S.V = double(S.V);
                S.c = double(S.c);
            else
                error("Only single or double precision arrays allowed. GpuArray/dlarray are coming.")
            end
        end

    end


    methods % get and check methods

        % get a box bounding the zonotope
        function B = getBox(obj)
            
            lb = cast(zeros(obj.dim, 1), 'like', obj.V);
            ub = cast(zeros(obj.dim, 1), 'like', obj.V);
                        
            for i=1:obj.dim
                lb(i) = obj.c(i) - norm(obj.V(i, :), 1);
                ub(i) = obj.c(i) + norm(obj.V(i, :), 1);
            end
            
            B = Box(lb, ub);
            
        end
        
        % return possible max indexes
        function max_ids = getMaxIndexes(obj)
            % @max_ids: index of the state that can be a max point
            
            % author: Dung Tran
            % date: 4/1/2020
            
            new_rs  = obj.toStar; 
            new_rs = new_rs.toImageStar(obj.dim, 1, 1);
            max_id = new_rs.get_localMax_index([1 1], [obj.dim 1], 1);
            max_ids = max_id(:, 1);

        end
        
        % check if a zonotope contains a point
        function bool = contains(obj, x)
            % @x: point
            % @bool: = '1' if the zonotope contain x else = '0'
            
            % author: Dung Tran
            % date: 1/9/2020
            
            if size(x,1) ~= obj.dim
                error('In consistent dimensions between the input point and the zonotope');
            end
            if size(x, 2) ~= 1
                error('Invalid input point, should have 1 column');
            end
            
            d = x - obj.c; 
            abs_V = abs(obj.V);
            d1 = sum(abs_V, 2); 
            
            x1 = (d <= d1);
            x2 = (d >= -d1);
            
            bool = (sum(x1) == obj.dim) && (sum(x2) == obj.dim);
                        
        end
        
        % get bounds of a zonotope
        function [lb, ub] = getBounds(obj)
            % clip method from StanleyBak
            % to get bound of a zonotope
            
            % Author: Dung Tran
            % Date: 1/4/2020
            
            pos_mat = obj.V'; 
            neg_mat = obj.V';
            pos_mat(pos_mat < 0) = 0; 
            neg_mat(neg_mat > 0) = 0;
            pos1_mat = cast(ones(1, size(obj.V, 2)), 'like', obj.V); 
            ub = transpose(pos1_mat*(pos_mat - neg_mat));
            lb = -ub;
            
            ub = obj.c + ub;
            lb = obj.c + lb;
            
        end
        
        % get ranges of a zonotope
        function [lb, ub] = getRanges(obj)
            B = obj.getBox;
            lb = B.lb;
            ub = B.ub;
        end
        
        % get range of a zonotope at specific index
        function [lb, ub] = getRange(obj, index)
            % @index: index of the state x[index] 
            % @lb: lower bound of x[index]
            % @ub: upper bound of x[index]
            
            % author: Dung Tran
            % date: 5/3/2019
            
            if index <= 0 || index > obj.dim
                error('Invalid index');
            end
            
            lb = obj.c(index) - norm(obj.V(index, :), 1);
            ub = obj.c(index) + norm(obj.V(index, :), 1);
        end
        
        % check if a index is larger than other
        function bool = is_p1_larger_than_p2(obj, p1_id, p2_id)
            % @p1_id: index of point 1
            % @p2_id: index of point 2
            % @bool = 1 if there exists the case that p1 >= p2
            %       = 0 if there is no case that p1 >= p2
            
            % author: Dung Tran
            % date: 7/10/2020
            
            S = obj.toStar;
            bool = S.is_p1_larger_than_p2(p1_id, p2_id);
        end
        
        % get an oriented rectangular hull enclosing a zonotope
        function Z = getOrientedBox(obj)
            % this code is from MATTISE of Prof. Girard, in 2005. 
            % date: 10/30/2018
            
            [Q, ~, ~] = svd(obj.V);
            P = Q' * obj.V;
            D = cast(diag(sum(abs(P), 2)), 'like', obj.V);
            Z = Zono(obj.c, Q*D);
            
        end
        
        % get interval hull
        function I = getIntervalHull(obj)
                      
            B = obj.getBox();
            I = B.toZono();
        end
        
        % get sup_{x\in Z }||x||_\infinity
        function r = getSupInfinityNorm(obj)
            
            B = obj.getBox();
            V1 = [abs(B.lb) abs(B.ub)];
            r = max(max(V1));
            
        end
        
        % get all vertices of a zonotope
        function V = getVertices(obj)
            
            % author: Dung Tran
            % date: 10/25/2018
            
            n = size(obj.V, 2); % number of generator
            
            lb = cast(-ones(n, 1), 'like', obj.V);
            ub = cast(ones(n, 1), 'like', obj.V);
            B = Box(lb, ub);
            V1 = B.getVertices();
            m = size(V1, 2); % number of vertices of the zonotope
            V = [];
            for i=1:m
                v = obj.c + obj.V * V1(:, i);
                V = [V v];
            end
            
        end
        
    end
    

    methods(Static) % plotting methods
        
        % plot an array of zonotopes
        function plots(varargin)
            % @Z: an array of zonotope
            % author: Dung Tran
            % date: 10/9/2019
            % update:4/2/2020
            
            switch nargin
                case 2
                    Z = varargin{1};
                    color = varargin{2};
                case 1
                    Z = varargin{1};
                    color = 'red';
                otherwise
                    error('Invalid number of input arguments');
            end
            
            n = length(Z);
            for i=1:n
                if ~isa(Z, 'Zono')
                    error('Z(%d) is not a zonotope', i);
                end
            end
            
            max_order = 12;
            
            for i=1:n-1
                max_gens = max_order * Z(i).dim;
                if max_gens >= size(Z(i).V, 2)
                    Zono.plot(Z(i), color);
                else
                    Zr = Z(i).orderReduction_box(max_gens);
                    Zono.plot(Zr,color);
                end
                hold on;
            end
            
            max_gens = max_order * Z(n).dim;
            if max_gens >= size(Z(n).V, 2)
                Zono.plot(Z(n),color);
            else
                Zr = Z(n).orderReduction_box(max_gens);
                Zono.plot(Zr,color);
            end
            
        end
        
        % plot a zonotope
        function plot(varargin)
            % author: Dung Tran
            % date: 10/20/2019
            % update: 4/2/2020
            
            switch nargin
                case 2
                    Z = varargin{1};
                    color = varargin{2};
                case 1
                    Z = varargin{1};
                    color = 'red';
                otherwise
                    error('Invalid number of input arguments');
            end
  
            P = Z.toPolyhedron();
            P.plot('color', color);
            
        end
    end
    
end

