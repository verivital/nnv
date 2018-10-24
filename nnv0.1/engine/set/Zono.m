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
        
        % convert to polyhedron
        function P = toPolyhedron(obj)
            
            n = size(obj.V, 2);           
            lb = -ones(n, 1);
            ub = ones(n, 1);                    
            Pa = Polyhedron('lb', lb, 'ub', ub);
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
        
    end
    
end

