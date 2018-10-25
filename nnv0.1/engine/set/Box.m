classdef Box
    %Hyper-rectangle class
    % Box and simple methods
    % Dung Tran: 10/4/2018
    properties
        lb = [];
        ub = [];
        center = []; 
        generators = []; 
    end
    
    methods
        
        % constructor
        function obj = Box(lb, ub)
            % @lb: lower-bound vector
            % @ub: upper-bound vector
            
            [n1, m1] = size(lb);
            [n2, m2] = size(ub);
            
            if m1 ~= 1 || m2 ~= 1
                error('lb and ub should be a vector');
            end
            
            if n1 ~= n2
                error('Inconsistent dimensions between lb and ub');
            end
            
            obj.lb = lb;
            obj.ub = ub;
            
            obj.center = 0.5 * (lb + ub);
            vec = 0.5 * (ub - lb);
            for i=1:n1
                if vec(i) ~= 0
                    gen = zeros(n1, 1);
                    gen(i) = vec(i);
                    obj.generators = [obj.generators gen];
                end                
            end
                
        end
        
        % affine mapping of a box
        
        function B = affineMap(obj, W, b)
            % @W: mapping matrix
            % @b: mapping vector
            
            % @B: return a new box bound the affine map
            
            if size(b, 2) ~= 1
                error('b should be a vector');
            end
            
            if size(W, 1) ~= size(b, 1) 
                error('Inconsistency between mapping matrix, mapping vector');
            end
            
            new_center = W * obj.center + b;
            new_generators = W * obj.generators;
            
            n = length(new_center);
            new_lb = zeros(n, 1);
            new_ub = zeros(n, 1);
            
            for i=1:n
                v = new_generators(i, :)';
                new_lb(i) = new_center(i) - norm(v, 1);
                new_ub(i) = new_center(i) + norm(v, 1);
            end
            
            
            B = Box(new_lb, new_ub);
            
        end
        
        
        % transform box to polyhedron
        function P = toPolyhedron(obj)
            
            P = Polyhedron('lb', obj.lb, 'ub', obj.ub);
        end
        
        % plot a box using mpt toolbox
        function plot(obj)
            P = obj.toPolyhedron;
            P.plot;
        end
        
        % transform box to star set
        function S = toStar(obj)
            
            c = (obj.lb + obj.ub)/2;
            V = obj.generators;
            n = size(obj.lb, 1);
            lb1 = -ones(n, 1);
            ub1 = ones(n, 1);
            P = Polyhedron('lb', lb1, 'ub', ub1);
            
            V = horzcat(c, V);
            S = Star(V, P.A, P.b);       
        end
        
        % transform box to zonotope
        function Z = toZono(obj)
            Z = Zono(obj.center, obj.generators);
        end
        
        % get all vertices of the box
        function V = getVertices(obj)
            % author: Dung Tran
            % date: 10/25/2018
            
            n = length(obj.lb);            
            N = 2^n; % number of vertices in the worst case 
            V = [];
            for i=0:N-1               
                b = de2bi(i, n+3);
                v = zeros(n, 1);
                for j=1:n
                    if b(j) == 1
                        v(j) = obj.ub(j);
                    else
                        v(j) = obj.lb(j);
                    end
                end                
                V = [V v];                
            end
            
            % delete duplicate vertices
            V1 = V';
            V1 = unique(V1, 'rows', 'stable');
            V = V1';
                      
        end
        
        
    end
    
    methods(Static)
       
        % box merging
        function B = boxHull(I)
            % merge boxes into one box
            % @I: array of boxes
            
            n = length(I);
            lb = [];
            ub = [];
            
            for i=1:n
                lb = [lb I.lb];
                ub = [ub I.ub];                
            end
            
            lb = min(lb, [], 2);
            ub = max(ub, [], 2);
            
            
            B = Box(lb, ub);
            
        end
        
        % plot an array of boxes 
        function plots(boxes)
            
            n = length(boxes);
            
            R = [];
            for i=1:n
                R = [R boxes(i).toPolyhedron];
            end
            R.plot;
            
        end
        
    end
    
end

