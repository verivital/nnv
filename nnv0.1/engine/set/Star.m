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
            
        end
        
        % affine mapping of star set S = Wx + b;
        function S = affineMap(obj, W, b)
            % @W: mapping matrix
            % @b: mapping vector
            % @S: new star set
            
            if size(W, 2) ~= obj.dim
                error('Inconsistency between the affine mapping matrix and dimension of the star set');
            end
            
            if size(b, 1) ~= size(W, 1)
                error('Inconsistency between the mapping vec and mapping matrix');
            end
            
            if size(b, 2) ~= 1
                error('Mapping vector should have one column');
            end
            
            newV = W * obj.V;
            newV(:, 1) = newV(:, 1) + b;
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
           
            V3 = horzcat(V1, V2);
            new_c = obj.V(:, 1) + X.V(:, 1);
            new_V = horzcat(new_c, V3);        
            new_C = blkdiag(obj.C, X.C);        
            new_d = vertcat(obj.d, X.d);
            
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
        
    end
    
end

