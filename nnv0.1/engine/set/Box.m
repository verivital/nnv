classdef Box
    %Hyper-rectangle class
    % Box and simple methods
    % Dung Tran: 10/4/2018
    properties
        lb = [];
        ub = [];
        dim = 0;
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
            obj.dim = length(lb);
            
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
                        
            c = 0.5 * (obj.lb + obj.ub);
            vec = 0.5 * (obj.ub - obj.lb);  
           
            V = [];
            alp_min = [];
            alp_max = [];
            for i=1:obj.dim
                if vec(i) ~= 0
                    gen = zeros(obj.dim, 1);
                    gen(i) = vec(i);
                    V = [V gen]; 
                    alp_min = [alp_min obj.lb(i)];
                    alp_max = [alp_max obj.ub(i)];
                                       
                end                
            end
            
            alp_min = alp_min';
            alp_max = alp_max';
           
            n = length(alp_min);
            C = vertcat(eye(n), -eye(n)); % constraint matrix
            d = vertcat(alp_max, -alp_min); % constraint vector
           
            V = horzcat(c, V);
            S = Star(V, C, d);       
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
        
         % plot list of boxes in 2D
        function plotBoxes_2D(boxes, x_pos, y_pos, color)
            % plot two dimension boxes
            
            
            n = length(boxes);
            
            for i=1:n
                if ~isa(boxes(i), 'Box')
                    error('plot object is not a box');
                end

                if x_pos > boxes(i).dim || y_pos > boxes(i).dim
                    error('Invalid x_pos or y_pos');
                end
                
                x = [boxes(i).lb(x_pos) boxes(i).ub(x_pos) boxes(i).ub(x_pos) boxes(i).lb(x_pos)];
                y = [boxes(i).lb(y_pos) boxes(i).lb(y_pos) boxes(i).ub(y_pos) boxes(i).ub(y_pos)];
                
                hold on;
                patch(x, y, color);
                                
            end
             
        end
        
        % plot list of boxes in 2D with no fill
        function plotBoxes_2D_noFill(boxes, x_pos, y_pos, color)
            n = length(boxes);
            
            for i=1:n
                if ~isa(boxes(i), 'Box')
                    error('plot object is not a box');
                end

                if x_pos > boxes(i).dim || y_pos > boxes(i).dim
                    error('Invalid x_pos or y_pos');
                end
                
                x = [boxes(i).lb(x_pos) boxes(i).ub(x_pos) boxes(i).ub(x_pos) boxes(i).lb(x_pos) boxes(i).lb(x_pos)];
                y = [boxes(i).lb(y_pos) boxes(i).lb(y_pos) boxes(i).ub(y_pos) boxes(i).ub(y_pos) boxes(i).lb(y_pos)];
                
                hold on;
                plot(x, y, color);
                                
            end
        end
        
        
        % plot list of boxes in 3D
        function plotBoxes_3D(boxes, x_pos, y_pos, z_pos, color)
            
            n = length(boxes);
            
            for i=1:n
                if ~isa(boxes(i), 'Box')
                    error('plot object is not a box');
                end

                if x_pos > boxes(i).dim || y_pos > boxes(i).dim || z_pos > boxes(i).dim
                    error('Invalid x_pos, y_pos or z_pos');
                end
                
                p1 = [boxes(i).lb(x_pos) boxes(i).lb(y_pos) boxes(i).lb(z_pos)];
                p2 = [boxes(i).ub(x_pos) boxes(i).lb(y_pos) boxes(i).lb(z_pos)];
                p3 = [boxes(i).ub(x_pos) boxes(i).lb(y_pos) boxes(i).ub(z_pos)];
                p4 = [boxes(i).lb(x_pos) boxes(i).lb(y_pos) boxes(i).ub(z_pos)];
                p5 = [boxes(i).ub(x_pos) boxes(i).ub(y_pos) boxes(i).lb(z_pos)];
                p6 = [boxes(i).ub(x_pos) boxes(i).ub(y_pos) boxes(i).ub(z_pos)];
                p7 = [boxes(i).lb(x_pos) boxes(i).ub(y_pos) boxes(i).ub(z_pos)];
                p8 = [boxes(i).lb(x_pos) boxes(i).ub(y_pos) boxes(i).lb(z_pos)];
                
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
        
        
    end
    
end

