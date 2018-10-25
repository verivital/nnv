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
    
    methods
        
        %constructor
        function obj = Zono(c, V)
            % @c: center vector
            % @V: generator matrix
            
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
        % generally not a zonotope, this function return an
        % over-approximation (a zonotope) of the convex hull.
        % main reference: 1) Reachability of Uncertain Linear Systems Using
        % Zonotopes, Antoin Girard, HSCC2005.
        function Z = convexHull(obj, X)
            % @X: input zonotope
            % @Z: output zonotope
            
            % author: Dung Tran
            % date: 10/25/2018
            
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
            %                  = (Z1 + e*Z1)/2 + (Z2 - e*Z2)/2  
            %                  where, '+' denote minkowski sum of two zonotopes
            % From minkowski sum method, one can see that:
            %    (Z1 + Z2)/2 =  0.5 * (c1 + c2, <g1, ..., gp, h1, ..., hq>)
            %    (Z1 - Z2)/2 = 0.5 * (c1 - c2, <g1, ..., gp, -h1, ...., -hq>)
            %    ea*(Z1-Z2)/2 = 0.5 * (0, <c1 - c2, g1, ..., gp, -h1, ..., -hq>)
            % Therefore:
            %     CH(Z1 U Z2) = 0.5(c1 + c2, <g1, ..., gp, h1, ..., hq, c1 - c2,
            %                                               g1, ..., gp, -h1, ..., -hq>)
            % So, the zonotope that over-approximate the convex hull of two zonotop
            % has 2(p + q) + 1 generators.
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
            new_V = 0.5 * [obj.V X.V (obj.c - X.c) obj.V -X.V];
            Z = Zono(new_c, new_V);
            
        end
        
             
        
        % convert to polyhedron
        function P = toPolyhedron(obj)
            
            n = size(obj.V, 2);           
            lb = -ones(n, 1);
            ub = ones(n, 1);                    
            B = Box(lb, ub);
            Vs = B.getVertices(); % get all vertices of a box
            Pa = Polyhedron('V', Vs');
            P = Pa.affineMap(obj.V, 'vrep') + obj.c;
        end
        
        % convert to Star
        function S = toStar(obj)
            n = size(obj.V, 2);           
            lb = -ones(n, 1);
            ub = ones(n, 1);                    
            Pa = Polyhedron('lb', lb, 'ub', ub);
            
            S = Star([obj.c obj.V], Pa.A, Pa.b);            
        end
        
        % plot a zonotope
        function plot(obj)
            P = obj.toPolyhedron();
            P.plot;
        end
        
        % get interval hull
        function B = getBox(obj)
            
            lb = zeros(obj.dim, 1);
            ub = zeros(obj.dim, 1);
            
            for i=1:obj.dim
                lb(i) = obj.c(i) - norm(obj.V(i, :), 1);
                ub(i) = obj.c(i) + norm(obj.V(i, :), 1);
            end
            
            B = Box(lb, ub);
            
        end
        
        % get all vertices of a zonotope
        function V = getVertices(obj)
            
            % author: Dung Tran
            % date: 10/25/2018
            
            n = size(obj.V, 2); % number of generator
            
            lb = -ones(n, 1);
            ub = ones(n, 1);
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
    
    methods(Static)
        
        % plot an array of zonotopes
        function plots(Z)
            % @Z: an array of zonotope
            
            n = length(Z);
            for i=1:n
                if ~isa(Z, 'Zono')
                    error('Z(%d) is not a zonotope', i);
                end
            end
            
            for i=1:n-1                
                Z(i).plot;
                hold on;
            end
            Z(n).plot;
            
        end
    end
    
end

