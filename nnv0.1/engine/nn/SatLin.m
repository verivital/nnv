classdef SatLin
    % SATLIN : class for computing reachable set of Satlin Transfer Function 
    %   Reference: https://www.mathworks.com/help/deeplearning/ref/satlin.html
    % Author: Dung Tran
    % Date: 27/2/2019
    
    properties
        
    end
    
    methods(Static) % evaluate method and reachability analysis with stars    
        
        % evaluation
        function y = evaluate(x)
            y = satlin(x);
        end
            
        % stepReach method, compute reachable set for a single step
        function S = stepReach(I, index)
            % @I: single star set input
            % @index: index of the neural performing stepSatLin
            % @S: star output set
            
            % author: Dung Tran
            % date: 27/2/2019
            
            
            if ~isa(I, 'Star')
                error('Input is not a star set');
            end
            
            % case 1: x(index) <= 0, SatLin(x(index)) = 0
            C0 = I.V(index, 2:I.nVar + 1);
            d0 = I.V(index, 1);
            C1 = [I.C; C0];
            d1 = [I.d; -d0];
            V1 = I.V;
            V1(index, :) = zeros(1, I.nVar + 1);
            S1 = Star(V1, C1, d1);
            if S1.isEmptySet
                S1 = [];
            end
            
            % case 2: 0 <= x(index) <= 1, SatLin(x(index)) = x(index)            
            V2 = I.V;
            C2 = [I.C; C0; -C0];
            d2 = [I.d; 1-d0; d0];
            S2 = Star(V2, C2, d2);
            if S2.isEmptySet
                S2 = [];
            end
            
            % case 3: x(index) >= 1
            C3 = [I.C; -C0];
            d3 = [I.d; d0 - 1];
            V3 = I.V;
            V3(index, 1) = 1;
            V3(index, 2:I.nVar + 1) = zeros(1, I.nVar);
            S3 = Star(V3, C3, d3);
            if S3.isEmptySet
                S3 = [];
            end
            
            S = [S1 S2 S3];
                     
        end
        
        
        % stepReach with multiple inputs
        function S = stepReachMultipleInputs(varargin)
            % @I: an array of stars
            % @index: index where stepReach is performed
            % @option: = 'parallel' use parallel computing
            %          = not declare -> don't use parallel computing
            
            % author: Dung Tran
            % date: 27/2/2019
            
            switch nargin
                case 3
                    I = varargin{1};
                    index = varargin{2};
                    option = varargin{3};
                case 2
                    I = varargin{1};
                    index = varargin{2};
                    option = [];
                otherwise
                    error('Invalid number of input arguments (should be 2 or 3)');
            end
            
            
            
            p = length(I);
            S = [];
            
            if isempty(option)
                
                for i=1:p
                    S =[S, SatLin.stepReach(I(i), index)];
                end
                
            elseif strcmp(option, 'parallel')
                
                parfor i=1:p
                    S =[S, SatLin.stepReach(I(i), index)];
                end
                
            else
                error('Unknown option');
            end
            
            
        end
       
        
        
        % function reachability analysis using Star
        function S = reach(varargin)
            % @I: an array of star input sets
            % @option: = 'parallel' use parallel option
            %          = '' do use parallel option
            
            % author: Dung Tran
            % date: 27/2/2019
            
            switch nargin
                case 2
                    I = varargin{1};
                    option = varargin{2};
                case 1
                    I = varargin{1};
                    option = [];
                otherwise
                    error('Invalid number of input arguments (should be 1 or 2)');
            end
            
            if ~isempty(I)       
                dim = I(1).dim;
                In = I;
                for i=1:dim
                    fprintf('\nPerforming SatLin_%d operation', i);
                    In = SatLin.stepReachMultipleInputs(In, i, option);
                end             
                
                S = In;
            else
                S = [];
            end
            
              
        end
        
        
        % step over-approximate reachability analysis using Star
        function S = stepReachStarApprox(I, index)
            % @I: star input set
            % @index: index of the neuron where we perform step Reach
            % @S: Star output set 
            
            % author: Dung Tran
            % date: 4/3/2019
            
            if ~isa(I, 'Star')
                error('Input is not a star');
            end
            
            [lb, ub] = I.getRange(index);
            
            if ub <= 0
                V = I.V;
                V(index, :) = zeros(1, I.nVar + 1);
                S = Star(V, I.C, I.d);
            end
            if lb >= 1
                V = I.V;
                V(index, :) = zeros(1, I.nVar + 1);
                V(index, 1) = 1;
                S = Star(V, I.C, I.d);                
            end
            if (1 > lb) && (lb > 0) && ub > 1
                % constraint 1: y(index) <= x[index]
                C1 = [-I.V(index, 2:I.nVar + 1) 1];
                d1 = I.V(index, 1);
                % constraint 2: y[index] <= 1
                C2 = zeros(1, I.nVar + 1);
                C2(I.nVar + 1) = 1;
                d2 = 1;
                % constraint 3: y[index] >= (1-lb)x/(ub-lb) + lb*(ub-1)/(ub-lb)
                C3 = [((1-lb)/(ub-lb))*I.V(index, 2:I.nVar+1) -1];
                d3 = -lb*(ub-1)/(ub-lb) - (1-lb)*I.V(index,1)/(ub-lb);
                
                m = size(I.C, 1);
                C0 = [I.C zeros(m, 1)];
                d0 = I.d;
                new_C = [C0; C1; C2; C3];
                new_d = [d0; d1; d2; d3];
                new_V = [I.V zeros(I.dim, 1)];
                new_V(index, :) = zeros(1, I.nVar+2);
                new_V(index, I.nVar+2) = 1;

                S = Star(new_V, new_C, new_d);
            end
            
            if lb >= 0 && ub <= 1
                S = Star(I.V, I.C, I.d);
            end
            if lb < 0 && (0 < ub) && (ub <= 1)
                
                n = I.nVar + 1;
                % over-approximation constraints 
                % constraint 1: y[index] = >= 0
                C1 = zeros(1, n);
                C1(n) = -1; 
                d1 = 0;
                % constraint 2: y[index] >= x[index]
                C2 = [I.V(index, 2:n) -1];
                d2 = -I.V(index, 1);
                % constraint 3: y[index] <= ub(x[index] -lb)/(ub - lb)
                C3 = [-(ub/(ub-lb))*I.V(index, 2:n) 1];
                d3 = -(ub*lb/(ub-lb))*(1 - I.V(index, 1));

                m = size(I.C, 1);
                C0 = [I.C zeros(m, 1)];
                d0 = I.d;
                new_C = [C0; C1; C2; C3];
                new_d = [d0; d1; d2; d3];
                new_V = [I.V zeros(I.dim, 1)];
                new_V(index, :) = zeros(1, n+1);
                new_V(index, n+1) = 1;

                S = Star(new_V, new_C, new_d);               
            end
            
            if lb < 0 && ub > 1
                % constraint 1: y[index] >= 0
                n = I.nVar + 1;
                C1 = zeros(1, n);
                C1(n) = -1; 
                d1 = 0;
                % constraint 2: y[index] <= 1
                C2 = zeros(1, n);
                C2(n) = 1;
                d2 = 1;
                % constraint 3: y[index] <= x/(1 -lb) - lb/(1-lb)
                C3 = [(-1/(1-lb))*I.V(index, 2:n) 1];
                d3 = (1/(1-lb))*I.V(index, 1) - lb/(1-lb);
                % constraint 4: y[index] >=  x/ub
                C4 = [(1/ub)*I.V(index, 2:n) -1];
                d4 = -(1/ub)*I.V(index, 1);
                
                m = size(I.C, 1);
                C0 = [I.C zeros(m, 1)];
                d0 = I.d;
                new_C = [C0; C1; C2; C3; C4];
                new_d = [d0; d1; d2; d3; d4];
                new_V = [I.V zeros(I.dim, 1)];
                new_V(index, :) = zeros(1, n+1);
                new_V(index, n+1) = 1;

                S = Star(new_V, new_C, new_d); 
                
            end
                 
        end
        
        % over-approximate reachability analysis using star
        function S = reach_star_approx(I)
            % @I: star input set
            % @S: star output set
            
            % author: Dung Tran
            % date: 4/3/2019
            
            if ~isa(I, 'Star')
                error('Input is not a star');
            end

            if isEmptySet(I)
                S = [];
            else
                In = I;
                for i=1:I.dim
                    fprintf('\nPerforming SatLin_%d operation', i);
                    In = SatLin.stepReachStarApprox(In, i);
                end
                S = In;
            end
        end
        
             
    end
    
    
    methods(Static) % over-approximate reachability analysis using zonotope
        
        % step reachability analysis using zonotope
        function Z = stepReachZonoApprox(I, index)
            % @I: zonotope input set
            % @index: index of the neuron performing step reach
            % @Z: zonotope output set
            
            % author: Dung Tran
            % date: 5/3/2019
            
            if ~isa(I, 'Zono')
                error('Input set is not a Zonotope');
            end
            
            [lb, ub] = I.getRange(index);
              
            if ub <= 0
                V = I.V;
                V(index, :) = zeros(1, size(V, 2));
                c = I.c;
                c(index) = 0;
                Z = Zono(c, V);
            end
            if lb >= 1
                V = I.V;
                V(index, :) = zeros(1, size(V,2));
                c = I.c;
                c(index) = 1;
                Z = Zono(c, V);                
            end
            
            if lb >= 0 && ub <= 1
                Z = Zono(I.c, I.V);
            end
            
            if (1 > lb) && (lb > 0) && ub > 1
                c = I.c;
                c(index) = c(index) + (1-ub)/2;
                V = zeros(I.dim, 1);
                V(index) = (1-ub)/2;
                V = [I.V V];
                Z = Zono(c, V);
            end
            
            if lb < 0 && (0 < ub) && (ub <= 1)
                
                c = I.c;
                lamda_opt = ub/(ub - lb);
                mu = -0.5*ub*lb/(ub - lb);
                c(index) = lamda_opt * c(index) + mu;
                V = zeros(I.dim, 1);
                V(index) = mu;
                V = [I.V V];
                Z = Zono(c, V);
                     
            end
            
            if lb < 0 && ub > 1
                
                % x + 1 -ub <= y <= x - lb, lb <= x <= ub
                c = I.c;
                mu = (1 + lb - ub)/2;
                c(index) = c(index) - lb  + mu;
                V = zeros(I.dim, 1);
                V(index) = mu;
                V = [I.V V];
                Z = Zono(c, V);
                
            end
                      
        end
        
        % over-approximate reachability analysis use zonotope
        function Z = reach_zono_approx(I)
            % @I: zonotope input
            % @Z: zonotope output
            
            % author: Dung Tran
            % date: 5/3/2019
            
            % reference: Fast and Effective Robustness Ceritification,
            % Gagandeep Singh, NIPS 2018
            
            if ~isa(I, 'Zono')
                error('Input is not a Zonotope');
            end
                      
            In = I;
            for i=1:I.dim
                fprintf('\nPerforming SatLin_%d operation', i);
                In = SatLin.stepReachZonoApprox(In, i);
            end
            Z = In;
            
        end
        
        
    end
    
    
    methods(Static) % reachability analysis using abstract-domain
        
        % future supporting method
        
    end
    
    methods(Static) % reachability analysis method using face-latice
        
        % future supporting method
        
    end
    
end

