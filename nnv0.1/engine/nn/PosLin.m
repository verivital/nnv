classdef PosLin
    % POSLIN Class contains method for reachability analysis for Layer with ReLU activation function 
    % IT GONNA BE USED TO REPLACE ReLU class in the future
    % Its name is consistent with the name of ReLU acivation function in
    % matlab called 'poslin'
    
    % author: Dung Tran
    % date: 27/2/2019
    
    properties
        
    end
    
    methods(Static) % evaluate method and reachability analysis with stars 
        
        % evaluation
        function y = evaluate(x)
            y = poslin(x);
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
            
            % case 2: 0 <= x(index), SatLin(x(index)) = x(index)            
            V2 = I.V;
            C2 = [I.C; -C0];
            d2 = [I.d; d0];
            S2 = Star(V2, C2, d2);
            if S2.isEmptySet
                S2 = [];
            end
                        
            S = [S1 S2];
                     
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
                    S =[S, PosLin.stepReach(I(i), index)];
                end
                
            elseif strcmp(option, 'parallel')
                
                parfor i=1:p
                    S =[S, PosLin.stepReach(I(i), index)];
                end
                
            else
                error('Unknown option');
            end
            
            
        end
        
        
        % exact reachability analysis using star
        function S = reach_star_exact(I, option)
            % @I: star input sets
            % @option: = 'parallel' using parallel computing
            %          = ''    do not use parallel computing
            
            % author: Dung Tran
            % date: 3/16/2019
            
             if ~isempty(I)
                dim = I(1).dim;
                In = I;
                for i=1:dim
                    fprintf('\nPerforming exact PosLin_%d operation using Star', i);
                    In = PosLin.stepReachMultipleInputs(In, i, option);
                end             
                S = In;
            else
                S = [];
            end
        end
        
              
        % step reach approximation using star
        function S = stepReachStarApprox(I, index)
            % @I: star set input
            % @index: index of the neuron performing stepReach
            % @S: star output

            % author: Dung Tran
            % date: 4/3/2019

            if ~isa(I, 'Star')
                error('Input is not a star');
            end

            [lb, ub] = I.getRange(index);

            if lb >= 0
                S = Star(I.V, I.C, I.d);
            elseif ub <= 0
                V = I.V;
                V(index, :) = zeros(1, I.nVar + 1);
                S = Star(V, I.C, I.d);
            elseif lb < 0 && ub > 0
                n = I.nVar + 1;
                % over-approximation constraints 
                % constraint 1: y[index] = ReLU(x[index]) >= 0
                C1 = zeros(1, n);
                C1(n) = -1; 
                d1 = 0;
                % constraint 2: y[index] >= x[index]
                C2 = [I.V(index, 2:n) -1];
                d2 = -I.V(index, 1);
                % constraint 3: y[index] <= ub(x[index] -lb)/(ub - lb)
                C3 = [-(ub/(ub-lb))*I.V(index, 2:n) 1];
                d3 = -ub*lb/(ub-lb) +  ub*I.V(index, 1)/(ub-lb);

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

        end


        % over-approximate reachability analysis using Star
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
                    fprintf('\nPerforming approximate PosLin_%d operation using Star', i);
                    In = PosLin.stepReachStarApprox(In, i);
                end
                S = In;
            end

        end



    end


    methods(Static) % over-approximate reachability analysis use zonotope
        
        % step over-approximate reachability analysis using zonotope
        function Z = stepReachZonoApprox(I, index)
            % @I: zonotope input set
            % @index: index of neuron performing stepReach
            % @Z: zonotope output set
            
            % author: Dung Tran
            % date: 5/3/2019
        
            % reference: Fast and Effective Robustness Ceritification,
            % Gagandeep Singh, NIPS 2018
           
            if ~isa(I, 'Zono')
                error('Input is not a Zonotope');
            end
               
            [lb, ub] = I.getRange(index);
            
            if lb >= 0
                Z = Zono(I.c, I.V);
                
            elseif ub <= 0
                c = I.c;
                c(index) = 0;
                V = I.V;
                V(index, :) = zeros(1, size(I.V, 2));
                Z = Zono(c, V);
                
            elseif lb < 0 && ub > 0
                
                lamda = ub/(ub -lb);
                mu = -0.5*ub*lb/(ub - lb);               
                
                c = I.c; 
                c(index) = lamda * c(index) + mu;
                V = I.V;
                V(index, :) = lamda * V(index, :);
                I1 = zeros(I.dim,1);
                I1(index) = mu;
                V = [V I1];
                
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
                fprintf('\nPerforming approximate PosLin_%d operation using Zonotope', i);
                In = PosLin.stepReachZonoApprox(In, i);
            end
            Z = In;
                       
        end
        
        
    end
    
    
    methods(Static) % reachability analysis using abstract-domain
        
        % step over-approximate reachability analysis using abstract-domain
        % we use star set to represent abstract-domain
        function S = stepReachAbstractDomain(I, index)
            % @I: star-input set
            % @index: index of neuron performing stepReach
            % @S: star output set represent abstract-domain of the output
            % set
            
            % author: Dung Tran
            % date: 5/3/2019
        
            % reference: An Abstract Domain for Certifying Neural Networks,
            % Gagandeep Singh, POPL 2019
           
            if ~isa(I, 'Star')
                error('Input is not a Star');
            end
               
            [lb, ub] = I.getRange(index);
            
            if lb >= 0
                S = Star(I.V, I.C, I.d);
            elseif ub <= 0
                V = I.V;
                V(index, :) = zeros(1, I.nVar + 1);
                S = Star(V, I.C, I.d);
            elseif lb < 0 && ub > 0
                
                S1 = ub*(ub-lb)/2; % area of the first candidate abstract-domain
                S2 = -lb*(ub-lb)/2; % area of the second candidate abstract-domain
                
                n = I.nVar + 1;
                if S1 < S2
                    % choose the first cadidate as the abstract-domain                  
                    % over-approximation constraints 
                    % constraint 1: y[index] = ReLU(x[index]) >= 0
                    C1 = zeros(1, n);
                    C1(n) = -1; 
                    d1 = 0;
                    % constraint 2: y[index] <= ub(x[index] -lb)/(ub - lb)
                    C2 = [-(ub/(ub-lb))*I.V(index, 2:n) 1];
                    d2 = -(ub*lb/(ub-lb)) + ub*I.V(index, 1)/(ub-lb);
                    
                    m = size(I.C, 1);
                    C0 = [I.C zeros(m, 1)];
                    d0 = I.d;
                    new_C = [C0; C1; C2];
                    new_d = [d0; d1; d2];
                    new_V = [I.V zeros(I.dim, 1)];
                    new_V(index, :) = zeros(1, n+1);
                    new_V(index, n+1) = 1;

                    S = Star(new_V, new_C, new_d);
                    
                else
                    % choose the second candidate as the abstract-domain                   
                    % over-approximation constraints 
                    % constraint 1: y[index] = ReLU(x[index]) >= x[index]
                    C1 = [I.V(index, 2:n) -1];
                    d1 = -I.V(index, 1);
                    % constraint 2: y[index] <= ub(x[index] -lb)/(ub - lb)
                    C2 = [-(ub/(ub-lb))*I.V(index, 2:n) 1];
                    d2 = -(ub*lb/(ub-lb)) + ub*I.V(index, 1)/(ub-lb);
                    m = size(I.C, 1);
                    C0 = [I.C zeros(m, 1)];
                    d0 = I.d;
                    new_C = [C0; C1; C2];
                    new_d = [d0; d1; d2];
                    new_V = [I.V zeros(I.dim, 1)];
                    new_V(index, :) = zeros(1, n+1);
                    new_V(index, n+1) = 1;

                    S = Star(new_V, new_C, new_d);
                                      
                end
                
            end
                       
        end
        
        
        % over-approximate reachability analysis using abstract-domain
        function S = reach_abstract_domain(I)
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
                    fprintf('\nPerforming approximate PosLin_%d operation using abstract domain', i);
                    In = PosLin.stepReachAbstractDomain(In, i);
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
            % date: 27/2/2019
            
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
                
                R = PosLin.reach_star_exact(I, option);
                
            elseif strcmp(method, 'approx-star')  % over-approximate analysis using star
                
                R = PosLin.reach_star_approx(I);
                
            elseif strcmp(method, 'approx-zono')  % over-approximate analysis using zonotope
                
                R = PosLin.reach_zono_approx(I);
                
            elseif strcmp(method, 'abs-dom')  % over-approximate analysis using abstract-domain
                
                R = PosLin.reach_abstract_domain(I);
                
            elseif strcmp(method, 'exact-face-latice') % exact analysis using face-latice
                fprintf('\nNNV have not yet support Exact Face-Latice Method');
            elseif strcmp(method, 'approx-face-latice') % over-approximate analysis using face-latice
                fprintf('\nNNV have not yet support Approximate Face-Latice Method');
            end
            
            
            
            
              
        end
        
        
    end
    
    
end

