classdef Star
    %Star set class 
    %   Star set defined by x = c + a[1]*v[1] + a[2]*v[2] + ... + a[n]*v[n]
    %                         = V * b, V = [c v[1] v[2] ... v[n]]
    %                                  b = [1 a[1] a[2] ... a[n]]^T                                   
    %                       where C*a <= d, constraints on a[i]
    %   Dung Tran: 10/22/2018
    
    properties
        V = []; % basic matrix
        C = []; % constraint matrix
        d = []; % constraint vector
        dim = 0; % dimension of star set
        nVar = 0; % number of variable in the constraints
        
        predicate_lb = []; % lower bound vector of predicate variable
        predicate_ub = []; % upper bound vector of predicate variable
                
    end
    
    methods
        
        % constructor
        function obj = Star(varargin)
            % @V: bassic matrix
            % @C: constraint matrix
            % @d: constraint vector
            
            switch nargin
                
                case 5
                    
                    V = varargin{1};
                    C = varargin{2};
                    d = varargin{3};
                    pred_lb = varargin{4};
                    pred_ub = varargin{5};
                    
                    [nV, mV] = size(V);
                    [nC, mC] = size(C);
                    [nd, md] = size(d);
                    [n1, m1] = size(pred_lb);
                    [n2, m2] = size(pred_ub);
                    
                    if mV ~= mC + 1
                        error('Inconsistency between basic matrix and constraint matrix');
                    end

                    if nC ~= nd
                        error('Inconsistency between constraint matrix and constraint vector');
                    end

                    if md ~= 1
                        error('constraint vector should have one column');
                    end
                    
                    if m1 ~=1 || m2 ~=1 
                        error('predicate lower- or upper-bounds vector should have one column');
                    end
                    
                    if n1 ~= n2 || n1 ~= mC
                        error('Inconsistency between number of predicate variables and predicate lower- or upper-bounds vector');
                    end

                    obj.V = V;
                    obj.C = C; 
                    obj.d = d;
                    obj.dim = nV;
                    obj.nVar = mC;
                    obj.predicate_lb = pred_lb;
                    obj.predicate_ub = pred_ub;
                
                case 3
                    V = varargin{1};
                    C = varargin{2};
                    d = varargin{3};
                    [nV, mV] = size(V);
                    [nC, mC] = size(C);
                    [nd, md] = size(d);

                    if mV ~= mC + 1
                        error('Inconsistency between basic matrix and constraint matrix');
                    end

                    if nC ~= nd
                        error('Inconsistency between constraint matrix and constraint vector');
                    end

                    if md ~= 1
                        error('constraint vector should have one column');
                    end

                    obj.V = V;
                    obj.C = C; 
                    obj.d = d;
                    obj.dim = nV;
                    obj.nVar = mC;
                                                            
                case 2
                    
                    % construct star from lower bound and upper bound
                    % vector
                    lb = varargin{1};
                    ub = varargin{2};
                    
                    B = Box(lb,ub);
                    S = B.toStar;
                    obj.V = S.V;
                    obj.C = S.C;
                    obj.d = S.d;
                    obj.dim = S.dim;
                    obj.nVar = S.nVar;
                    obj.predicate_lb = -ones(S.nVar, 1);
                    obj.predicate_ub = ones(S.nVar, 1);
                    
                
                case 0
                    % create empty Star (for preallocation an array of star)
                    obj.V = [];
                    obj.C = [];
                    obj.d = [];
                    obj.dim = 0;
                    obj.nVar = 0;
                    
                otherwise
                    
                    error('Invalid number of input arguments (should be 0 or 2 or 3 or 5)');
            end
            
            
        end
        
        % check is empty set
        function bool = isEmptySet(obj)
            
            options = optimoptions(@linprog, 'Preprocess', 'none', 'Display','none');           
            f = ones(1, obj.nVar);
            
            [~, ~, exitflag, ~] = linprog(f, obj.C, obj.d, [], [], [], [], [], options);
            
            if exitflag == 1
                bool = 0;
            else
                bool = 1;
            end 
            
        end
        
        % check if a star set is a subset of other
        function bool = isSubSet(obj, S)
            
            if ~isa(S, 'Star')
                error('Input set is not a Star');
            end
            
            if obj.dim ~= S.dim
                error('Two sets have different dimension');
            end
            
            P1 = obj.toPolyhedron();
            P2 = S.toPolyhedron();
            
            bool = (P1 <= P2);
            
        end
        
        % check if a star set contain a point
        function bool = contains(obj, s)
            % @s: a star point
            % @bool: = 1 star set contains s, else no
            
            % author: Dung Tran
            % date: 1/3/2019
            
            
            if size(s,1) ~= obj.dim
                error('Dimension mismatch');
            end
            if size(s,2) ~= 1
                error('Invalid star point');
            end
            
            A = obj.C;
            b = obj.d;
            Ae = obj.V(:, 2:obj.nVar+1);
            be = s - obj.V(:,1);
            
            P = Polyhedron('A', A, 'b', b, 'Ae', Ae, 'be', be);
            
            bool = ~P.isEmptySet;
                     
        end
        
        % sampling a star set
        function V = sample(obj, N)
            % @N: number of points in the samples
            % @V: a set of at most N sampled points in the star set 
            
            % author: Dung Tran
            % date: 1/3/2019
            
            if N < 1
                error('Invalid number of samples');
            end
            
            B = obj.getBox;
            if isempty(B)
                V = [];
            else
                lb = B.lb;
                ub = B.ub;
                
                X = cell(1, obj.dim);
                V1 = [];
                for i=1:obj.dim
                    X{1, i} = (ub(i) - lb(i)).*rand(2*N, 1) + lb(i);
                    V1 = vertcat(V1, X{1, i}');
                end
                
                V = [];
                for i=1:2*N
                    if obj.contains(V1(:, i))
                        V = [V V1(:, i)];
                    end
                end
                
                if size(V, 2) > N
                    V = V(:, 1:N);
                end             
                
            end
                     
        end
        
        % affine mapping of star set S = Wx + b;
        function S = affineMap(obj, W, b)
            % @W: mapping matrix
            % @b: mapping vector
            % @S: new star set
            
            if size(W, 2) ~= obj.dim
                error('Inconsistency between the affine mapping matrix and dimension of the star set');
            end
            
            if ~isempty(b)
                if size(b, 1) ~= size(W, 1)
                    error('Inconsistency between the mapping vec and mapping matrix');
                end

                if size(b, 2) ~= 1
                    error('Mapping vector should have one column');
                end

                newV = W * obj.V;
                newV(:, 1) = newV(:, 1) + b;
            else
                newV = W * obj.V;
            end
            
            if ~isempty(obj.predicate_lb)
                S = Star(newV, obj.C, obj.d, obj.predicate_lb, obj.predicate_ub);
            else
                S = Star(newV, obj.C, obj.d);
            end
           
        end
        
        % Minkowski Sum
        function S = MinkowskiSum(obj, X)
            % @X: another star with same dimension
            % @S: new star
            
            if ~isa(X, 'Star')
                error('Input matrix is not a Star');
            else
                if X.dim ~= obj.dim
                    error('Input star and current star have different dimensions');
                end
            end
            
            m1 = size(obj.V, 2);
            m2 = size(X.V, 2);
            
            V1 = obj.V(:, 2:m1);
            V2 = X.V(:, 2:m2);
            new_c = obj.V(:, 1) + X.V(:, 1);
            
            % check if two Star has the same constraints
            if size(obj.C, 1) == size(X.C, 1) && size(obj.C, 2) == size(X.C, 2) && norm(obj.C - X.C) + norm(obj.d - X.d) < 0.0001
                  
               V3 = V1 + V2;
               new_V = horzcat(new_c, V3);
               S = Star(new_V, obj.C, obj.d); % new Star has the same number of basic vectors
             
            else
                
                V3 = horzcat(V1, V2);
                new_c = obj.V(:, 1) + X.V(:, 1);
                new_V = horzcat(new_c, V3);        
                new_C = blkdiag(obj.C, X.C);        
                new_d = vertcat(obj.d, X.d);

                S = Star(new_V, new_C, new_d); % new Star has more number of basic vectors
                
            end
                
        end
        
        % intersection with a half space: H(x) := Hx <= g
        function S = intersectHalfSpace(obj, H, g)
            % @H: HalfSpace matrix
            % @g: HalfSpace vector
            % @S : new star set with more constraints
            
            % author: Dung Tran
            % date: 10/27/2018
            
            [nH, mH] = size(H);
            [ng, mg] = size(g);
            
            if mg ~= 1
                error('Halfspace vector should have one column');
            end
            if nH ~= ng
                error('Inconsistent dimension between Halfspace matrix and halfspace vector');
            end
            if mH ~= obj.dim
                error('Inconsistent dimension between halfspace and star set');
            end
            
            m = size(obj.V, 2);
            
            C1 = H*obj.V(:, 2:m);
            d1 = g - H*obj.V(:,1);
            
            new_C = vertcat(obj.C, C1);
            new_d = vertcat(obj.d, d1);
                       
            S = Star(obj.V, new_C, new_d);
            
            if S.isEmptySet
                S = [];
            end
            
        end
        
        % scalar map of a Star S' = alp * S, 0 <= alp <= alp_max
        function S = scalarMap(obj, alp_max)
            % @a_max: maximum value of a
            % @S: new Star
            
            % note: we always require that alp >= 0
            % author: Dung Tran
            % date: 10/27/2018
            
            % =============================================================
            % S: x = alp*c + V* alph * a, Ca <= d
            % note that:   Ca <= d -> C*alph*a <= alp*a <= alp_max * d
            % let: beta = alp * a, we have
            % S := x = alp * c + V * beta, C * beta <= alp_max * d,
            %                              0 <= alp <= alp_max
            % Let g = [beta; alp]
            % S = Star(new_V, new_C, new_d), where:
            %   new_V = [0 c V], new_C = [0 -1; 0 1; 0 C], new_d = [0; alpha_max; alp_max * d]
            %       
            % S has one more basic vector compared with obj
            % =============================================================
            
            new_c = zeros(obj.dim, 1);
            new_V = [obj.V new_c];
            new_C = blkdiag(obj.C, [-1; 1]);           
            new_d = vertcat(alp_max*obj.d, 0, alp_max); 
            S = Star(new_V, new_C, new_d);           
        end
        
        
        % convex hull of two Stars
        function S = convexHull(obj, X)
            % @X: input star
            % @S: an over-approximation of (convex hull) of two Stars
            
            % author: Dung Tran
            % date: 10/27/2018
            
            % =============================================================
            % This method is similar to the method proposed by Prof. Girard
            % for Zonotope. See the following paper:
            % 1) Reachability of Uncertain Linear Systems Using
            % Zonotopes, Antoin Girard, HSCC2005.
            %
            % CH(S1 U S2) := {l*x1 + (1-l)*x2}, 0 <= l <= 1, x1 \in S1, x2
            % \in S2.     := {x2 + l*(x1 - x2)}
            % 
            %             := c2 + a2 * V2 + l(c1 - c2) + l*a1*V1 - l*a2*V2
            % let beta1 = l * a1, beta2 = l*a2 -> C1 * beta1 <= l * d1 <=
            % d1 and C2 * beta2 <= l * d2 <= d2
            %
            % CH(S1 U S2) = c2 + a2 * V2 + l * (c1 - c2) + beta1 * V1 -
            %                                                   beta2 * V2
            %
            % if S2 = L*S1, i.e., S2 is a projection of S1
            % 
            % CH(S1 U S2) := {(l + (1-l)*L) * x1} = {(l(I - L) + L) * x1}
            %             := (l * (I - L) * x1 + L * x1)
            % =============================================================
            
            if ~isa(X, 'Star')
                error('Input set is not a Star');
            end
            
            if X.dim ~= obj.dim
                error('Inconsisten dimension between input set and this star');
            end

            %new_V = horzcat(X.V(:,1), X.V(:, 2:X.nVar+1), (obj.V(:,1) - X.V(:,1)), obj.V(:,2:obj.nVar+1), -X.V(:, 2:X.nVar+1));
            %new_C = blkdiag(X.C, [-1; 1], obj.C, X.C);
            %new_d = vertcat(X.d, [0; 1], obj.d, X.d);      
            
            %S = Star(new_V, new_C, new_d);
            
            c1 = obj.V(:, 1);
            c2 = X.V(:, 1);
            V1 = obj.V(:, 2:obj.nVar + 1);
            V2 = X.V(:, 2:obj.nVar + 1);
            
            new_V = horzcat(c1-c2, V1, -V2);
            new_C = blkdiag(obj.C, X.C);
            new_d = vertcat(obj.d, X.d);
            
            S2 = Star(new_V, new_C, new_d);
            S = X.MinkowskiSum(S2);
                                    
        end
        
        % convexHull of a Star with its linear transformation
        function S = convexHull_with_linearTransform(obj, L)
            % @L: linear transformation matrix
            % @S: new Star
            
            % author: Dung Tran
            % date: 10/27/2018
            
            [nL, mL] = size(L);
            
            if nL ~= mL 
                error('Trasformation matrix should be a square matrix');
            end
            if nL ~= obj.dim
                error('Inconsistent dimension of tranformation matrix and this zonotope');
            end
            
            M = eye(obj.dim) - L; 
            S1 = obj.affineMap(L, []);
            S2 = obj.affineMap(M, []);
            
            new_V = [S1.V(:,1) S1.V(:,2:S1.nVar+1) S2.V(:,1) S2.V(:,2:S2.nVar+1)];
            new_C = blkdiag(S1.C, [-1; 1], S2.C);
            new_d = vertcat(S1.d, [0; 1], S2.d);
            
            S = Star(new_V, new_C, new_d);
                        
        end
        
        % order reduction for Stars
        % similar for order reduction for zonotope
        % see the Zono class for the detail
        % We reduce the number of basic vectors of a Star
        function S = orderReduction_box(obj, n_max)
            % @n_max: maximum allowable number of basic vectors
            % @S: a new Star with number of basic vectors = n_max
            
            % author: Dung Tran - not finish yet
            % date: 10/29/2018
            
            if n_max < obj.dim
                error('n_max should be >= %d', obj.dim);
            end
            
            if n_max == obj.dim
                B = obj.getBox();
                S = B.toStar();
            end
            
            n = size(obj.V, 2) - 1; % number of generators
            if n_max > obj.dim
                n1 = n - n_max + obj.dim; % number of generators need to be reduced
                                
                % sort generators based on their lengths
                length_gens = zeros(n, 1);
                for i=2:n+1
                    length_gens(i) = norm(obj.V(:, i), 2);
                end
                [~, sorted_ind] = sort(length_gens); 
                
                sorted_gens = zeros(obj.dim, n);
                for i=1:n
                    sorted_gens(:, i) = obj.V(:, sorted_ind(i));
                end
                
                S = [];
                
            end
            
        end
        
        
        % convert to polyhedron
        function P = toPolyhedron(obj)
            
            b = obj.V(:, 1);        
            W = obj.V(:, 2:size(obj.V, 2));           
            Pa = Polyhedron('A', obj.C, 'b', obj.d);
            P = Pa.affineMap(W) + b;
            
        end
        
        % convert to 2D Star set
        function S = toStar2D(obj, height, width)
            % @height: height of 2D Star set
            % @width: width of 2D Star set
            % @S: a 2D Star set, please see Star2D class
            
            % author: Dung Tran
            % date: 12/20/2018
            
            if height * width ~= obj.dim
                error('Height and Width of 2D Star should satisfy height * width = dimension of current 1D star set');
            end
            
            new_V = cell(1, obj.nVar + 1);
            for i=1:obj.nVar + 1
                new_V{i} = transpose(reshape(obj.V(:,i), [height, width]));
            end
            
            S = Star2D(new_V, obj.C, obj.d);
            
        end
        
        % plot star set
        function plot(obj)
            P = obj.toPolyhedron();
            P.plot;
        end
        
        % find a box bounding a star
        function B = getBox(obj)
            
            lb = zeros(obj.dim, 1);
            ub = zeros(obj.dim, 1); 

            if isempty(obj.C) || isempty(obj.d) % star set is just a vector (one point)
                lb = obj.V(:, 1);
                ub = obj.V(:, 1);

            else % star set is a set

                options = optimset('Display','none');

                for i=1:obj.dim
                    f = obj.V(i, 2:obj.nVar + 1);
                    %[~, fval, exitflag, ~] = linprog(f, obj.C, obj.d, [], [], [], [], [], options);
                    [~, fval, exitflag, ~] = glpk(f, obj.C, obj.d, [], [], [], [], [], options);
                    if exitflag > 0
                        lb(i) = fval + obj.V(i, 1);
                    else
                        lb = [];
                        ub = [];
                        break;
                    end
                    %[~, fval, exitflag, ~] = linprog(-f, obj.C, obj.d, [], [], [], [], [], options);
                    [~, fval, exitflag, ~] = glpk(-f, obj.C, obj.d, [], [], [], [], [], options);
                    if exitflag > 0
                        ub(i) = -fval + obj.V(i, 1);
                    else
                        lb = [];
                        ub = [];
                        break;
                    end

                end

            end         
            
            if isempty(lb) || isempty(ub)
                B = [];
            else
                B = Box(lb, ub);           
            end
            
        end
        
        % find a box bounding a star
        function B = getBox_parallel(obj)          

            if isempty(obj.C) || isempty(obj.d) % star set is just a vector (one point)
                lb = obj.V(:, 1);
                ub = obj.V(:, 1);

            else % star set is a set
                lb = zeros(obj.dim,1);
                ub = zeros(obj.dim,1); 
                options = optimset('Display','none');
                parfor i=1:obj.dim
                   
                    f = obj.V(i, 2:obj.nVar + 1);
                    %[~, fval, exitflag, ~] = linprog(f, obj.C, obj.d, [], [], [], [], [], options);
                    [~, fval, exitflag, ~] = glpk(f, obj.C, obj.d, [], [], [], [], [], options);
                    if exitflag > 0        
                        lb(i) = fval + obj.V(i,1);
                    else
                        lb(i) = -Inf;                         
                    end
                    %[~, fval, exitflag, ~] = linprog(-f, obj.C, obj.d, [], [], [], [], [], options);
                    [~, fval, exitflag, ~] = glpk(-f, obj.C, obj.d, [], [], [], [], [], options);
                    if exitflag > 0
                        ub(i) = -fval + obj.V(i, 1);
                    else
                        ub(i) = Inf;
                    end

                end

            end         
            
            if isempty(lb) || isempty(ub)
                B = [];
            else
                B = Box(lb, ub);           
            end
            
        end
        
        % estimates ranges of state vector quickly
        % these ranges are not the exact range, it is an
        % over-approximation of the exact ranges
        function B = getBoxFast(obj)
            
            % author: Dung Tran
            % date: 5/29/2019
            
            if isempty(obj.C) || isempty(obj.d) % star set is just a vector (one point)
                lb = obj.V(:, 1);
                ub = obj.V(:, 1);
                B = Box(lb, ub);
            else      
                [pred_lb, pred_ub] = obj.getPredicateBounds;         
                B1 = Box(pred_lb, pred_ub);
                Z = B1.toZono;
                Z = Z.affineMap(obj.V(:,2:obj.nVar + 1), obj.V(:,1));
                B = Z.getBox;
            end
                         
        end
        
        % get bounds of predicate variables
        function [pred_lb, pred_ub] = getPredicateBounds(obj)
            % author: Dung Tran
            % date: 5/30/2019
            
            if isempty(obj.predicate_lb) || isempty(obj.predicate_ub)
                P = Polyhedron('A', obj.C, 'b', obj.d);
                P.outerApprox;
                pred_lb = P.Internal.lb;
                pred_ub = P.Internal.ub;
            else
                pred_lb = obj.predicate_lb;
                pred_ub = obj.predicate_ub;
            end  
            
        end
        
        
        % find range of a state at specific position
        function [xmin, xmax] = getRange(obj, index)
            % @index: position of the state
            % range: min and max values of x[index]
            
            % author: Dung Tran
            % date: 11/16/2018
            
            if index < 1 || index > obj.dim
                error('Invalid index');
            end
            
            options = optimset('Display','none');

            f = obj.V(index, 2:obj.nVar + 1);

            %[~, fval, exitflag, ~] = linprog(f, obj.C, obj.d, [], [], [], [], [], options);
            [~, fval, exitflag, ~] = glpk(f, obj.C, obj.d, [], [], [], [], [], options);
            if exitflag > 0
                xmin = fval + obj.V(index, 1);
            else
                error('Cannot find an optimal solution, exitflag = %d', exitflag);
            end          

            %[~, fval, exitflag, ~] = linprog(-f, obj.C, obj.d, [], [], [], [], [], options);
            [~, fval, exitflag, ~] = glpk(-f, obj.C, obj.d, [], [], [], [], [], options);
            if exitflag > 0
                xmax = -fval + obj.V(index, 1);
            else
                error('Cannot find an optimal solution');
            end
            
        end
        
        % find a oriented box bounding a star
        function B = getOrientedBox(obj)
            % author: Dung Tran
            % date: 11/8/2018
            
            [Q, Z, P] = svd(obj.V(:, 2:obj.nVar + 1));
            S = Z * P';
            lb = zeros(obj.dim, 1);
            ub = zeros(obj.dim, 1);
            
            options = optimset('Display','none');
            
            for i=1:obj.dim
                f = S(i, :);
                %[~,fval, exitflag, ~] = linprog(f, obj.C, obj.d, [], [], [], [], [], options);
                [~,fval, exitflag, ~] = glpk(f, obj.C, obj.d, [], [], [], [], [], options);
                
                if exitflag > 0
                    lb(i) = fval;
                else
                    error('Cannot find an optimal solution');
                end
                
                
                %[~, fval, exitflag, ~] = linprog(-f, obj.C, obj.d, [], [], [], [], [], options);
                [~,fval, exitflag, ~] = glpk(-f, obj.C, obj.d, [], [], [], [], [], options);
                if exitflag > 0
                    ub(i) = -fval;
                else
                    error('Cannot find an optimal solution');
                end
                
                
            end
            
            new_V = [obj.V(:,1) Q];
            new_C = vertcat(eye(obj.dim), -eye(obj.dim));
            new_d = vertcat(ub, -lb);
            
            B = Star(new_V, new_C, new_d);
            
        end
        
        % find a zonotope bounding a star (an over-approximation of a star using zonotope)
        function Z = getZono(obj)
            
            % author: Dung Tran
            % date: 10/25/2018
            
            B = obj.getBox;
            if ~isempty(B)
                Z = B.toZono;
            else
                Z = [];
            end
                        
        end
        
        % concatenate with other star
        function S = concatenate(obj, X)
            % @X: input star
            % @S: output star with is the concatenation of obj and X
            
            % author: Dung Tran
            % date: 11/18/2018
            
            if ~isa(X, 'Star')
                error('Input is not a Star');
            end
            
            c1 = obj.V(:, 1);
            c2 = X.V(:, 1);
            
            V1 = obj.V(:, 2:obj.nVar+1);
            V2 = X.V(:, 2:X.nVar + 1);
            
            c = vertcat(c1, c2);
            
            V1 = vertcat(V1, zeros(X.dim, obj.nVar));
            V2 = vertcat(zeros(obj.dim, X.nVar), V2);
            
            new_V = horzcat(c, V1, V2);
            new_C = blkdiag(obj.C, X.C);
            new_d = vertcat(obj.d, X.d);
            
            S = Star(new_V, new_C, new_d);            
            
        end
        
    end
    
    methods(Static) % plot methods
        
        % plot an array of Star (plot exactly, this is time consuming)
        function plots(S)
            % @S: an array of Stars
            
            n = length(S);
            P = [];
            for i=1:n
                P = [P S(i).toPolyhedron];
            end
            P.plot;
            
        end
        
        % plot ranges of stars at specific dimension versus time series
        function plotRanges_2D(stars, index, times, color)
            % @stars: an array of stars
            % @index: index of the specified dimension 
            % @times: an array of time points
            % @color: color of the plot
            
            % author: Dung Tran
            % date: 11/21/2018
            
            if length(stars) ~= length(times)
                error('Inconsistent length between stars and times');
            end
            
            n = length(stars);
            [ymin, ymax] = stars(1).getRange(index);
            y = [ymin ymax];
            x = [times(1) times(1)];
            plot(x, y, color);
            hold on;
            
            for i=2:n
                [ymin, ymax] = stars(i).getRange(index);
                y = [ymin ymin ymax ymax ymin];
                x = [times(i-1) times(i) times(i) times(i-1) times(i-1)];
                plot(x, y, color)
                hold on; 
            end
            hold off;
            
        end
        
        % plot an array of Stars using 2D boxes
        % plot list of boxes in 2D
        function plotBoxes_2D(stars, x_pos, y_pos, color)
            % plot stars using two dimension boxes
            % author: Dung Tran
            % date: 11/16/2018
                     
            n = length(stars);
            xmin = zeros(n,1);
            xmax = zeros(n,1);
            ymin = zeros(n,1);
            ymax = zeros(n,1);
            
            for i=1:n
                if isa(stars(i), 'Star')
                    [xmin(i), xmax(i)] = stars(i).getRange(x_pos);
                    [ymin(i), ymax(i)] = stars(i).getRange(y_pos);
                else
                    error('%d th input object is not a star', i);
                end
            end
            
                        
            for i=1:n
                                
                x = [xmin(i) xmax(i) xmax(i) xmin(i)];
                y = [ymin(i) ymin(i) ymax(i) ymax(i)];
                
                hold on;
                patch(x, y, color);
                                
            end
             
        end
        
        % plot list of boxes in 2D with no fill
        function plotBoxes_2D_noFill(stars, x_pos, y_pos, color)
            % plot stars using two dimension boxes without filling
            % author: Dung Tran
            % date: 11/16/2018
            
            n = length(stars);
            xmin = zeros(n,1);
            xmax = zeros(n,1);
            ymin = zeros(n,1);
            ymax = zeros(n,1);
            
            for i=1:n
                if isa(stars(i), 'Star')
                    [xmin(i), xmax(i)] = stars(i).getRange(x_pos);
                    [ymin(i), ymax(i)] = stars(i).getRange(y_pos);
                else
                    error('%d th input object is not a star', i);
                end
            end
            
                        
            for i=1:n
                                
                x = [xmin(i) xmax(i) xmax(i) xmin(i) xmin(i)];
                y = [ymin(i) ymin(i) ymax(i) ymax(i) ymin(i)];
                
                hold on;
                plot(x, y, color);
                                
            end
            
        end
        
        
        % plot list of boxes in 3D
        function plotBoxes_3D(stars, x_pos, y_pos, z_pos, color)
            % plot stars using three dimensionall boxes
            % author: Dung Tran
            % date: 11/16/2018
            
            
            n = length(stars);
            xmin = zeros(n,1);
            xmax = zeros(n,1);
            ymin = zeros(n,1);
            ymax = zeros(n,1);
            zmin = zeros(n,1);
            zmax = zeros(n,1);
            
            for i=1:n
                if isa(stars(i), 'Star')
                    [xmin(i), xmax(i)] = stars(i).getRange(x_pos);
                    [ymin(i), ymax(i)] = stars(i).getRange(y_pos);
                    [zmin(i), zmax(i)] = stars(i).getRange(z_pos);
                else
                    error('%d th input object is not a star', i);
                end
            end
            
            
            for i=1:n
                                
                p1 = [xmin(i) ymin(i) zmin(i)];
                p2 = [xmax(i) ymin(i) zmin(i)];
                p3 = [xmax(i) ymin(i) zmax(i)];
                p4 = [xmin(i) ymin(i) zmax(i)];
                p5 = [xmax(i) ymax(i) zmin(i)];
                p6 = [xmax(i) ymax(i) zmax(i)];
                p7 = [xmin(i) ymax(i) zmax(i)];
                p8 = [xmin(i) ymax(i) zmin(i)];
                
                % line p1->p2->p3 ->p4->p1
                                
                p = vertcat(p1, p2, p3, p4, p1);
                p = p';
                x = p(1,:);
                y = p(2, :);
                z = p(3, :);
                          
                               
                plot3(x, y, z, color);
                
                
                % line p5->p6->p7->p8->p5
               
                
                p = vertcat(p5, p6, p7, p8, p5);
                p = p';
                x = p(1,:);
                y = p(2, :);
                z = p(3, :);
                
                
                hold on;                
                plot3(x, y, z, color);
                
                % line p4->p7
                
                p = vertcat(p4, p7);
                p = p';
                x = p(1,:);
                y = p(2, :);
                z = p(3, :);
                
                
                hold on;                
                plot3(x, y, z, color);
                
                
                % line p3->p6
                p = vertcat(p3, p6);
                p = p';
                x = p(1,:);
                y = p(2, :);
                z = p(3, :);
                
                
                hold on;                
                plot3(x, y, z, color);
                
                % line p1->p8
                p = vertcat(p1, p8);
                p = p';
                x = p(1,:);
                y = p(2, :);
                z = p(3, :);

                hold on;                
                plot3(x, y, z, color);
                
                % line p2->p5
                
                p = vertcat(p2, p5);
                p = p';
                
                x = p(1,:);
                y = p(2, :);
                z = p(3, :);
                
               
                hold on;                
                plot3(x, y, z, color);
                
                hold on;
                                
            end
            
            hold off;
            
        end
        
        
        % bound a set of stars by a box
        function B = get_hypercube_hull(stars)
            % @stars: an array of stars
            % @S: a box (represented as a star)
            
            n = length(stars);
            
            dim = stars(1).dim;
            
            for i=2:n
                
                if ~isa(stars(i), 'Star')
                    error('The %d th object is not a star', i);
                end                
                if stars(i).dim ~= dim
                    error('Inconsistent dimensions between stars');
                end
                
            end
            
            lb = [];
            ub = [];
            for i=1:n
                B1 = stars(i).getBox();
                lb = [lb B1.lb];
                ub = [ub B1.ub];
            end
            
            lb = min(lb, [], 2);
            ub = max(ub, [], 2);
            
            B = Box(lb, ub);            
            
        end
        
        % convex hull of stars
        function S = get_convex_hull(stars)
            % @stars: an array of stars
            % @S: a new star which is a convex hull of the input stars
            
            % author: Dung Tran
            % date: 1/21/2019
            
            n = length(stars);
            dim = stars(1).dim; 
            
            
            for i=2:n
                if ~isa(stars(i), 'Star')
                    error('The %d th object is not a star', i);
                end
                if stars(i).dim ~= dim
                    error('Inconsistent dimensions between stars');
                end
            end
            
            for i=1:n
                X(i) = stars(i).toPolyhedron;
            end
            
            U = PolyUnion(X);
            S = U.convexHull;
            
        end
        
        % concatenate many stars
        
        function S = concatenateStars(stars)
            % @stars: an array of stars
            
            new_c = [];
            new_V = [];
            new_C = [];
            new_d = [];
            
            n = length(stars);
            
            for i=1:n
                if ~isa(stars(i), 'Star')
                    error('The %d th input is not a Star', i);
                end
                
                new_c = vertcat(new_c, stars(i).V(:,1));
                new_V = blkdiag(new_V, stars(i).V(:, 2:stars(i).nVar + 1));
                new_C = blkdiag(new_C, stars(i).C);
                new_d = vertcat(new_d, stars(i).d);
                
            end
            
            S = Star([new_c new_V], new_C, new_d);
           
        end
       
        
        % merge stars using boxes and overlapness
        function S = merge_stars(I, nS, parallel)
            % @I: array of stars
            % @nP: number of stars of the output S
            
            % author: Dung Tran
            % date: 2/25/2019
            
            n = length(I);
            B = [];
            if strcmp(parallel, 'single')             
                
                for i=1:n
                    B = [B I(i).getBox];
                end

                m = I(1).dim;

                n = length(B);

                C = zeros(n, 2*m);
                for i=1:n
                    C(i, :) = [B(i).lb' B(i).ub'];
                end

                idx = kmeans(C, nS); % clustering boxes into nP groups

                R = cell(nS, 1);

                for i=1:nS
                    for j=1:n
                        if idx(j) == i
                            R{i, 1} = [R{i, 1} B(j)];
                        end
                    end
                end

                S = [];
                for i=1:nS
                    B = Box.boxHull(R{i, 1}); % return a box                    
                    S = [S B.toStar];
                end
                
            elseif strcmp(parallel, 'parallel')
                
                parfor i=1:n
                    B = [B I(i).getBox];
                end

                m = I(1).dim;

                n = length(B);

                C = zeros(n, 2*m);
                for i=1:n
                    C(i, :) = [B(i).lb' B(i).ub'];
                end

                idx = kmeans(C, nS); % clustering boxes into nP groups

                R = cell(nS, 1);

                for i=1:nS
                    for j=1:n
                        if idx(j) == i
                            R{i, 1} = [R{i, 1} B(j)];
                        end
                    end
                end

                S = [];
                parfor i=1:nS
                    B = Box.boxHull(R{i, 1}); % return a box                    
                    S = [S B.toStar];
                end
                
            else
                error('Unknown parallel option');
            end

        end
           
    end
    
end

