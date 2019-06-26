classdef SatLins
    % SATLIN : class for computing reachable set of Satlins Transfer Function 
    %   Reference: https://www.mathworks.com/help/deeplearning/ref/satlins.html
    % Author: Dung Tran
    % Date: 27/2/2019
    
    properties
        
    end
    
    methods(Static) % evaluate method and reachability analysis with stars 
        
        % evaluation
        function y = evaluate(x)
            y = satlins(x);
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
            
            % case 1: x(index) <= -1, SatLin(x[index]) = -1
            C0 = I.V(index, 2:I.nVar + 1);
            d0 = I.V(index, 1);
            C1 = [I.C; C0];
            d1 = [I.d; -1-d0];
            V1 = I.V;
            V1(index, :) = zeros(1, I.nVar + 1);
            V1(index, 1) = -1;
            S1 = Star(V1, C1, d1);
            if S1.isEmptySet
                S1 = [];
            end
            
            % case 2: -1 <= x(index) <= 1, SatLin(x[index]) = x[index]            
            V2 = I.V;
            C2 = [I.C; C0; -C0];
            d2 = [I.d; 1-d0; d0+1];
            S2 = Star(V2, C2, d2);
            if S2.isEmptySet
                S2 = [];
            end
            
            % case 3: x(index) >= 1, SatLins(x[index]) = 1
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
                    S =[S, SatLins.stepReach(I(i), index)];
                end
                
            elseif strcmp(option, 'parallel')
                
                parfor i=1:p
                    S =[S, SatLins.stepReach(I(i), index)];
                end
                
            else
                error('Unknown option');
            end
            
            
        end
        
        
        % function reachability analysis using Star
        function S = reach_star_exact(I, option)
            % @I: an array of star input sets
            % @option: = 'parallel' use parallel option
            %          = '' do use parallel option
            
            % author: Dung Tran
            % date: 27/2/2019
            
            if ~isempty(I)       
                dim = I(1).dim;
                In = I;
                for i=1:dim
                    fprintf('\nPerforming SatLins_%d operation', i);
                    In = SatLins.stepReachMultipleInputs(In, i, option);
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
            
            if ub <= -1
                V = I.V;
                V(index, :) = zeros(1, I.nVar + 1);
                V(index, 1) = -1;
                S = Star(V, I.C, I.d);
            end
            
            if lb >= 1
                V = I.V;
                V(index, :) = zeros(1, I.nVar + 1);
                V(index, 1) = 1;
                S = Star(V, I.C, I.d);                
            end
            
            if lb >= -1 && ub <= 1
                S = Star(I.V, I.C, I.d);
            end
            
            
            if (1 > lb) && (lb > -1) && ub > 1
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
            

            if lb < -1 && (-1 < ub) && (ub <= 1)
                
                n = I.nVar + 1;
                % over-approximation constraints 
                % constraint 1: y[index] = >= -1
                C1 = zeros(1, n);
                C1(n) = -1; 
                d1 = 1;
                % constraint 2: y[index] >= x[index]
                C2 = [I.V(index, 2:n) -1];
                d2 = -I.V(index, 1);
                % constraint 3: y[index] <= ((ub+1)x[index] -ub(lb + 1))/(ub - lb)
                C3 = [-((ub + 1)/(ub -lb))*I.V(index, 2:n) 1];
                d3 = -ub*(lb + 1)/(ub - lb) + ((ub+1)/(ub-lb))*I.V(index,1);

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
            
            if lb < -1 && ub > 1
                % constraint 1: y[index] >= -1
                n = I.nVar + 1;
                C1 = zeros(1, n);
                C1(n) = -1; 
                d1 = 1;
                % constraint 2: y[index] <= 1
                C2 = zeros(1, n);
                C2(n) = 1;
                d2 = 1;
                % constraint 3: y[index] <= 2x/(1 -lb) - (lb+1)/(1-lb)
                C3 = [(-2/(1-lb))*I.V(index, 2:n) 1];
                d3 = (2/(1-lb))*I.V(index, 1) - (lb+1)/(1-lb);
                % constraint 4: y[index] >=  2x/(ub + 1) + (1-ub)/(1 + ub)
                C4 = [(2/(ub+1))*I.V(index, 2:n) -1];
                d4 = -(2/(ub+1))*I.V(index, 1) + (ub-1)/(ub+1);
                
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
                    fprintf('\nPerforming SatLins_%d operation', i);
                    In = SatLins.stepReachStarApprox(In, i);
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
              
            if ub <= -1
                V = I.V;
                V(index, :) = zeros(1, size(V, 2));
                c = I.c;
                c(index) = -1;
                Z = Zono(c, V);
            end
            if lb >= 1
                V = I.V;
                V(index, :) = zeros(1, size(V,2));
                c = I.c;
                c(index) = 1;
                Z = Zono(c, V);                
            end
            
            if lb >= -1 && ub <= 1
                Z = Zono(I.c, I.V);
            end
            
            if (1 > lb) && (lb >= -1) && ub > 1
                % x + 1 -ub <= y <= x, lb <= x <= ub
                c = I.c;
                c(index) = c(index) + (1-ub)/2;
                V = zeros(I.dim, 1);
                V(index) = (1-ub)/2;
                V = [I.V V];
                Z = Zono(c, V);
            end
            
            if lb < -1 && (0 < ub) && (ub <= 1)
                % lamda * x + lamda - 1 <= y <= lamda * x + ub*(1-lamda)
                % lamda_opt = ub/(ub - lb);
                c = I.c;
                lamda_opt = ub/(ub - lb);
                mu = (ub + 1)*(1 - lamda_opt)/2;
                c(index) = lamda_opt * c(index) + (lamda_opt - 1) + mu;
                V = zeros(I.dim, 1);
                V(index) = mu;
                V = [I.V V];
                Z = Zono(c, V);
                     
            end
            
            if lb < -1 && ub > 1
                % x + 1 -ub <= y <= x -1 - lb, lb <= x <= ub
                c = I.c;
                mu = (ub - lb - 2)/2;
                c(index) = c(index) + 1 - ub + mu;
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
            % date: 6/3/2019
            
            % reference: Fast and Effective Robustness Ceritification,
            % Gagandeep Singh, NIPS 2018
            
            if ~isa(I, 'Zono')
                error('Input is not a Zonotope');
            end
                      
            In = I;
            for i=1:I.dim
                fprintf('\nPerforming SatLins_%d operation', i);
                In = SatLins.stepReachZonoApprox(In, i);
            end
            Z = In;
            
        end

        
    end
    
    
    methods(Static) % reachability analysis using abstract-domain
        
        % step over-approximate reachability analysis using abstract-domain
        % we extend abstract-domain method from ReLU to SatLin
        % we use star set to represent abstract-domain
        function S = stepReachAbstractDomain(I, index)
            % @I: star-input set
            % @index: index of neuron performing stepReach
            % @S: star output set represent abstract-domain of the output
            % set
            
            % author: Dung Tran
            % date: 16/3/2019
        
            % reference: An Abstract Domain for Certifying Neural Networks,
            % Gagandeep Singh, POPL 2019
           
            if ~isa(I, 'Star')
                error('Input is not a Star');
            end
               
            [lb, ub] = I.getRange(index);
            
            
            if ub <= -1
                V = I.V;
                V(index, :) = zeros(1, I.nVar + 1);
                V(index, 1) = -1;
                S = Star(V, I.C, I.d);
            end
            
            if lb >= 1
                V = I.V;
                V(index, :) = zeros(1, I.nVar + 1);
                V(index, 1) = 1;
                S = Star(V, I.C, I.d);                
            end
            
            if lb >= -1 && ub <= 1
                S = Star(I.V, I.C, I.d);
            end
            
            
            if (1 > lb) && (lb > -1) && ub > 1
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
        
                new_V = [I.V zeros(I.dim, 1)];
                new_V(index, :) = zeros(1, I.nVar+2);
                new_V(index, I.nVar+2) = 1;
                
                S1 = (ub-lb)*(1-lb)/2; % area of the first candidate abstract domain
                S2 = (ub-lb)*(ub-1)/2; % area of the second candidate abstract domain
                
                if S1 < S2
                    % get first cadidate as resulted abstract-domain
                    new_C = [C0; C2; C3];
                    new_d = [d0; d2; d3];
                else
                    % get second candidate as resulted abstract-domain
                    new_C = [C0; C1; C3];
                    new_d = [d0; d1; d3];
                end
                
                S = Star(new_V, new_C, new_d);
            end
            

            if lb < -1 && (-1 < ub) && (ub <= 1)
                
                n = I.nVar + 1;
                % over-approximation constraints 
                % constraint 1: y[index] = >= -1
                C1 = zeros(1, n);
                C1(n) = -1; 
                d1 = 1;
                % constraint 2: y[index] >= x[index]
                C2 = [I.V(index, 2:n) -1];
                d2 = -I.V(index, 1);
                % constraint 3: y[index] <= ((ub+1)x[index] -ub(lb + 1))/(ub - lb)
                C3 = [-((ub + 1)/(ub -lb))*I.V(index, 2:n) 1];
                d3 = -ub*(lb + 1)/(ub - lb) + ((ub+1)/(ub-lb))*I.V(index,1);

                m = size(I.C, 1);
                C0 = [I.C zeros(m, 1)];
                d0 = I.d;

                new_V = [I.V zeros(I.dim, 1)];
                new_V(index, :) = zeros(1, n+1);
                new_V(index, n+1) = 1;
                
                S1 = (ub-lb)*(1+ub)/2; % area of the first candidate abstract domain
                S2 = (ub-lb)*(-1-lb)/2; % area of the second candidate abstract domain
                
                if S1 < S2
                    % get first cadidate as resulted abstract-domain
                    new_C = [C0; C1; C3];
                    new_d = [d0; d1; d3];
                else
                    % get second candidate as resulted abstract-domain
                    new_C = [C0; C2; C3];
                    new_d = [d0; d2; d3];
                end
                
                S = Star(new_V, new_C, new_d);               
            end
            
            if lb < -1 && ub > 1
                % constraint 1: y[index] >= -1
                n = I.nVar + 1;
                C1 = zeros(1, n);
                C1(n) = -1; 
                d1 = 1;
                % constraint 2: y[index] <= 1
                C2 = zeros(1, n);
                C2(n) = 1;
                d2 = 1;
                % constraint 3: y[index] <= 2x/(1 -lb) - (lb+1)/(1-lb)
                C3 = [(-2/(1-lb))*I.V(index, 2:n) 1];
                d3 = (2/(1-lb))*I.V(index, 1) - (lb+1)/(1-lb);
                % constraint 4: y[index] >=  2x/(ub + 1) + (1-ub)/(1 + ub)
                C4 = [(2/(ub+1))*I.V(index, 2:n) -1];
                d4 = -(2/(ub+1))*I.V(index, 1) + (ub-1)/(ub+1);
                
                m = size(I.C, 1);
                C0 = [I.C zeros(m, 1)];
                d0 = I.d;

                new_V = [I.V zeros(I.dim, 1)];
                new_V(index, :) = zeros(1, n+1);
                new_V(index, n+1) = 1;
                
                S1 = (ub-lb)*(1 + (2*ub-lb-1)/(1-lb))/2; % area of the first candidate abstract domain
                S2 = (ub-lb)*(1-((2*lb + 1 - ub)/(ub + 1)))/2; % area of the second candidate abstract domain
                
                if S1 < S2
                    % get first cadidate as resulted abstract-domain
                    new_C = [C0; C1; C3];
                    new_d = [d0; d1; d3];
                else
                    % get second candidate as resulted abstract-domain
                    new_C = [C0; C2; C4];
                    new_d = [d0; d2; d4];
                end
                
                S = Star(new_V, new_C, new_d); 
                
            end
            
            
                       
        end
        
        
        % over-approximate reachability analysis using abstract-domain
        function S = reach_abstract_domain(I)
            % @I: star input set
            % @S: star output set

            % author: Dung Tran
            % date: 16/3/2019


            if ~isa(I, 'Star')
                error('Input is not a star');
            end

            if isEmptySet(I)
                S = [];
            else
                In = I;
                for i=1:I.dim
                    fprintf('\nPerforming PosLin_%d operation', i);
                    In = SatLins.stepReachAbstractDomain(In, i);
                end
                S = In;
            end

        end     
        
        
    end
    
    methods(Static) % reachability analysis method using face-latice
        
        % future supporting method
        
    end
    
    
    methods(Static) % main reach method
        
        % main function for reachability analysis
        function R = reach(varargin)
            % @I: an array of star input sets
            % @method: 'exact-star' or 'approx-star' or 'approx-zono' or
            % 'abs-dom' or 'face-latice'
            % @option: = 'parallel' use parallel option
            %          = '' do use parallel option
            
            % author: Dung Tran
            % date: 16/3/2019
            
            switch nargin
                
                case 3
                    I = varargin{1};
                    method = varargin{2};
                    option = varargin{3};
                
                case 2
                    I = varargin{1};
                    method = varargin{2};
                    option = [];
                case 1
                    I = varargin{1};
                    method = 'exact-star';
                    option = [];
                otherwise
                    error('Invalid number of input arguments (should be 1, 2 or 3)');
            end
            
            
            if strcmp(method, 'exact-star') % exact analysis using star
                
                R = SatLins.reach_star_exact(I, option);
                
            elseif strcmp(method, 'approx-star')  % over-approximate analysis using star
                
                R = SatLins.reach_star_approx(I);
                
            elseif strcmp(method, 'approx-zono')  % over-approximate analysis using zonotope
                
                R = SatLins.reach_zono_approx(I);
                
            elseif strcmp(method, 'abs-dom')  % over-approximate analysis using abstract-domain
                
                R = SatLins.reach_abstract_domain(I);
                
            elseif strcmp(method, 'exact-face-latice') % exact analysis using face-latice
                fprintf('\nNNV have not yet support Exact Face-Latice Method');
            elseif strcmp(method, 'approx-face-latice') % over-approximate analysis using face-latice
                fprintf('\nNNV have not yet support Approximate Face-Latice Method');
            end
                        
        end
        
    end
    
end

