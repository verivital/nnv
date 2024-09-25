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
    
    methods % constructor and main methods
        
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

            % Creating a large matrix may result in out of memory errors using first method...
            try 
                % Speeding up implementation
                gens = diag(vec); % generate matrix
                if numel(gens) > 1
                    gens(:,all(gens==0)) = []; % delete colums with no info
                end
                obj.generators = gens;
            catch
                % This works well for large input sets with few perturbed pixels
                template_c = cast(zeros(obj.dim, 1), 'like', ub);
                gen_locs = find(vec~= 0);
                for i = gen_locs
                    gen = template_c;
                    gen(i) = vec(i);
                    obj.generators = [obj.generators gen];
                end

            end
            
        end
        
        % single partition of a box
        function Bs = singlePartition(obj, part_id, part_num)
            % @part_id: index of the state being partitioned
            % @part_num: number of partition
            
            % author: Dung Tran
            % date: 10/19/2019
            
            if part_id < 1 || part_id > obj.dim
                error('Invalid partition index');
            end
            
            if part_num < 1
                error('Invalid partition number');
            end
            
            if part_num == 1
                Bs = obj;
            else
                del = (obj.ub(part_id) - obj.lb(part_id))/part_num;
                
                Bs = [];
                
                for i=1:part_num
                    new_lb = obj.lb;
                    new_ub = obj.ub;                    
                    new_lb(part_id) = obj.lb(part_id) + (i-1)*del;
                    new_ub(part_id) = new_lb(part_id) + del;
                    Bs = [Bs Box(new_lb, new_ub)];
                end
                
            end

        end
        
        % partition a box into smaller boxes
        function Bs = partition(obj, part_indexes, part_numbers)
            % @part_indexes: the indexes of the states being
            %                    partitioned, an 1D array
            % @part_numbers: the number of partitions at specific
            %                    index, an 1D array
            
            % author: Dung Tran
            % date: 10/19/2019
            
            [n1, m1] = size(part_indexes); 
            [n2, m2] = size(part_numbers); 
            
            if (n1 ~= 1 || n2 ~= 1)
                error('Invalid part_indexes or part_numbers array, should have one row');
            end
            
            if m1 ~= m2
                error('Inconsistent between part_indexes and part_number array, should have the same number of elements');
            end
            
            for i=1:m1
                if ((part_indexes(i) < 1) || (part_indexes(i) > obj.dim))
                    error('The %d^th index in part_indexes array is invalid', i);
                end
                                
                if part_numbers(i) < 1
                    error('The %d^th number in part_numbers array is invalid', i);
                end                
            end
            
            Bs = [];
            B1 = obj; 
            for i=1:m1             
                n = length(B1);
                B2 = [];
                for j=1:n
                    B2 = [B2 B1(j).singlePartition(part_indexes(i), part_numbers(i))];
                end
                B1 = B2;                
            end
            
            Bs = B1;

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
            new_lb = cast(zeros(n, 1), 'like', W);
            new_ub = cast(zeros(n, 1), 'like', W);
            
            for i=1:n
                v = new_generators(i, :)';
                new_lb(i) = new_center(i) - norm(v, 1);
                new_ub(i) = new_center(i) + norm(v, 1);
            end
            
            B = Box(new_lb, new_ub);
            
        end

        % Generate samples from Box set
        function samples = sample(obj, N)
            % @N: number of points in the samples
            % @V: a set of at most N sampled points in the star set 
            
            % author: Diego Manzanas
            % date: 02/17/2023
            
            if N < 1
                error('Invalid number of samples');
            end
            n = obj.dim;
                        
            samples = (obj.ub - obj.lb).*rand(n,N) + obj.lb;
                
        end
        
    end


    methods % transformation methods

        % transform box to polyhedron
        function P = toPolyhedron(obj)
            P = Polyhedron('lb', obj.lb, 'ub', obj.ub);
        end

        % transform box to star set
        function S = toStar(obj)
            Z = obj.toZono;
            S = Z.toStar;        
        end
        
        % transform box to zonotope
        function Z = toZono(obj)
            Z = Zono(obj.center, obj.generators);
        end
    
    end


    methods % get methods

        % get Range
        function [lb, ub] = getRange(obj)
            lb = obj.lb;
            ub = obj.ub;
        end
        
        % get all vertices of the box
        function V = getVertices(obj)
            % author: Dung Tran
            % date: 10/25/2018
            
            n = length(obj.lb);            
            N = 2^n; % number of vertices in the worst case 
            V = [];
            for i=0:N-1
                b=dec2bin(i, n+3);
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

    end

    methods(Static) % plot methods
        
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
        
        % plot a box using mpt toolbox
        function plot(obj)
            P = obj.toPolyhedron;
            P.plot;
        end
    
    end
    
end

