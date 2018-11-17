classdef PLOT
    %PLOT Class contains static methods for efficiently plotting reachable set
    %   Dung Tran: 11/15/2018
    
    properties
        
    end
    
    methods(Static)
        
        % plot Polyhedra, Star or Boxes using boxes
        function plot_2d_box(sets, x_pos, y_pos, color)
           % @sets: array of input sets, can be Star, Boxes or Polyhedra
           % @x_pos: x position
           % @y_pos: y position
           % @color: color 
           
           % author: Dung Tran
           % date: 11/16/2018
           
           n = length(sets);
           xmin = zeros(n,1);
           xmax = zeros(n,1);
           ymin = zeros(n,1);
           ymax = zeros(n,1);
           
           for i=1:n
               
               if isa(sets(i), 'Polyhedron')
                   
                   [xmin(i), xmax(i)] = Conversion.getRange(sets(i), x_pos);
                   [ymin(i), ymax(i)] = Conversion.getRange(sets(i), y_pos);
                   
               elseif isa(sets(i), 'Star')
                   
                   [xmin(i), xmax(i)] = sets(i).getRange(x_pos);
                   [ymin(i), ymax(i)] = sets(i).getRange(y_pos);  
                   
               elseif isa(sets(i), 'Box')
                   
                   xmin(i) = sets(i).lb(x_pos);
                   xmax(i) = sets(i).ub(x_pos);
                   ymin(i) = sets(i).lb(y_pos);
                   ymax(i) = sets(i).ub(y_pos);
                   
               else
                   error('Unknown type of input set');
               end
               
           end
           
           boxes = [];
           for i=1:n
               lb = [xmin(i); ymin(i)];
               ub = [xmax(i); ymax(i)];
               boxes = [boxes Box(lb, ub)];
           end
           
           PLOT.plotBoxes_2D(boxes, 1, 2, color);
           
        end

        % plot Polyhedra, Star or Boxes using boxes
        function plot_2d_box_noFill(sets, x_pos, y_pos, color)
           % @sets: array of input sets, can be Star, Boxes or Polyhedra
           % @x_pos: x position
           % @y_pos: y position
           % @color: color 
           
           % author: Dung Tran
           % date: 11/16/2018
           
           n = length(sets);
           xmin = zeros(n,1);
           xmax = zeros(n,1);
           ymin = zeros(n,1);
           ymax = zeros(n,1);
           
           for i=1:n
               if isa(sets(i), 'Polyhedron')
                   
                   [xmin(i), xmax(i)] = Conversion.getRange(sets(i), x_pos);
                   [ymin(i), ymax(i)] = Conversion.getRange(sets(i), y_pos);
                   
               elseif isa(sets(i), 'Star')
                   [xmin(i), xmax(i)] = sets(i).getRange(x_pos);
                   [ymin(i), ymax(i)] = sets(i).getRange(y_pos);                  
               elseif isa(sets(i), 'Box')
                   
                   xmin(i) = sets(i).lb(x_pos);
                   xmax(i) = sets(i).ub(x_pos);
                   ymin(i) = sets(i).lb(y_pos);
                   ymax(i) = sets(i).ub(y_pos);             
                   
               else
                   error('Unknown type of input set');
               end
           end
           
           boxes = [];
           for i=1:n
               lb = [xmin(i); ymin(i)];
               ub = [xmax(i); ymax(i)];
               boxes = [boxes Box(lb, ub)];
           end
           
           PLOT.plotBoxes_2D_noFill(boxes, 1, 2, color);
           
        end
        
        % plot Polyhedra, Star or Boxes using 3d boxes
        function plot_3d_box(sets, x_pos, y_pos, z_pos, color)
            % @sets: array of input sets, can be Star, Boxes or Polyhedra
            % @x_pos: x position
            % @y_pos: y position
            % @z_pos: z position
            
            % @color: color 
            % author: Dung Tran
            % date: 11/16/2018
            
           n = length(sets);
           xmin = zeros(n,1);
           xmax = zeros(n,1);
           ymin = zeros(n,1);
           ymax = zeros(n,1);
           zmin = zeros(n,1);
           zmax = zeros(n,1);
           
           for i=1:n
               if isa(sets(i), 'Polyhedron')
                   
                   [xmin(i), xmax(i)] = Conversion.getRange(sets(i), x_pos);
                   [ymin(i), ymax(i)] = Conversion.getRange(sets(i), y_pos);
                   [zmin(i), zmax(i)] = Conversion.getRange(sets(i), z_pos);
                   
               elseif isa(sets(i), 'Star')
                   [xmin(i), xmax(i)] = sets(i).getRange(x_pos);
                   [ymin(i), ymax(i)] = sets(i).getRange(y_pos);
                   [zmin(i), zmax(i)] = sets(i).getRange(z_pos); 
                   
               elseif isa(sets(i), 'Box')
                   xmin(i) = sets(i).lb(x_pos);
                   xmax(i) = sets(i).ub(x_pos);
                   ymin(i) = sets(i).lb(y_pos);
                   ymax(i) = sets(i).ub(y_pos);
                   zmin(i) = sets(i).lb(z_pos);
                   zmax(i) = sets(i).ub(z_pos);
                   
               else
                   error('Unknown type of input set');
               end
           end
           
           boxes = [];
           for i=1:n
               lb = [xmin(i); ymin(i); zmin(i)];
               ub = [xmax(i); ymax(i); zmax(i)];
               boxes = [boxes Box(lb, ub)];
           end
           
           PLOT.plotBoxes_3D(boxes, 1, 2, 3, color);
          
            
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

