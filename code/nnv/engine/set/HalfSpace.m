classdef HalfSpace
    % HalfSpace class defining Gx <= g
    % Dung Tran: 3/27/2019
    
    properties
        G = []; % half-space matrix
        g = []; % half-space vector
        dim = 0; % dimension of half-space
    end
    
    methods
        
        % constructor
        function obj = HalfSpace(G, g)
            % @G: Half-space matrix 
            % @g: Half-space vector
           
            [n1, m1] = size(G);
            [n2, m2] = size(g);
            
            if n1 ~= n2
                error('Inconsistent dimension between half-space matrix and half-space vector');
            end
            
            if m2 ~= 1
                error('Half-space vector should have one column');
            end
            
            obj.G = G;
            obj.g = g;
            obj.dim = m1;
        end
        
        % check contain
        
        function bool = contains(obj, x)
            % @x: input vector
            % @bool: = 1 -> half-space contain point x
            %        = 0 -> half-space does not contain point x
            
            % author: Dung Tran
            % date: 3/27/2019
            
            [n, m] = size(x);
                        
            if n ~= obj.dim
                error('Inconsistent dimension between the vector x and the half-space object');
            end
            
            if m ~= 1
                error('Input vector x should have one column');
            end
            
            y = obj.g - obj.G * x;
            
            map = find(y < 0, 1); % which index has y(i) < 0
            
            bool = isempty(map);        
        
        end
        
        % plot half-space
        function plot(obj)
            % plot half-space
            
            P = Polyhedron('A', obj.G, 'b', obj.g);
            P.plot;
        end
        
    end
end

