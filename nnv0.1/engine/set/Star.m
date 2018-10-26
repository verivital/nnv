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
    end
    
    methods
        
        % constructor
        function obj = Star(V, C, d)
            % @V: bassic matrix
            % @C: constraint matrix
            % @d: constraint vector
            
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
            S = Star(newV, obj.C, obj.d);            
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

            new_V = horzcat(X.V(:,1), X.V(:, 2:X.nVar+1), (obj.V(:,1) - X.V(:,1)), obj.V(:,2:obj.nVar+1), -X.V(:, 2:X.nVar+1));
            display(new_V);
            new_C = blkdiag(X.C, [-1; 1], obj.C, X.C);
            display(new_C);
            new_d = vertcat(X.d, [0; 1], obj.d, X.d);
            
            S = Star(new_V, new_C, new_d);
                                    
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
        
        
        % convert to polyhedron
        function P = toPolyhedron(obj)
            
            b = obj.V(:, 1);        
            W = obj.V(:, 2:size(obj.V, 2));           
            Pa = Polyhedron('A', obj.C, 'b', obj.d);
            P = Pa.affineMap(W, 'vrep') + b;
            
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
            
            for i=1:obj.dim
                f = obj.V(i, 2:size(obj.V, 2));
                [~, fval] = linprog(f, obj.C, obj.d);
                lb(i) = fval + obj.V(i, 1);
                [~, fval] = linprog(-f, obj.C, obj.d);
                ub(i) = -fval + obj.V(i, 1);
            end           
            B = Box(lb, ub);           
        end
        
        % find a zonotope bounding a star (an over-approximation of a star using zonotope)
        function Z = getZono(obj)
            
            % author: Dung Tran
            % date: 10/25/2018
            
            P = Polyhedron('A', obj.C, 'b', obj.d);
            P.outerApprox;
            lb = P.Internal.lb;
            ub = P.Internal.ub;
            
            n = length(lb);
            new_V = [];
            new_c = obj.V(:,1);
            for i=1:n
                new_c = new_c + 0.5 * (ub(i) + lb(i)) * obj.V(:, i+1);
                new_V = [new_V 0.5 * (ub(i) - lb(i)) * obj.V(:, i+1)];
            end
            
            Z = Zono(new_c, new_V);
                        
        end
        
    end
    
    methods(Static)
        
        % plot an array of Star
        function plots(S)
            % @S: an array of Stars
            
            n = length(S);
            for i=1:n-1
                S(i).plot;
                hold on;
            end
            S(n).plot;
            
        end
        
    end
    
end

