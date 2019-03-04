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
                    fprintf('\nPerforming PosLin_%d operation', i);
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
                fprintf('\nPerforming PosLin_%d operation', i);
                In = PosLin.stepReachStarApprox(In, i);
            end
            S = In;
        end
              
    end
   
        
        
    end
    
    
    methods(Static) % reachability analysis using abstract-domain
        
        % future supporting method
        
    end
    
    methods(Static) % reachability analysis method using face-latice
        
        % future supporting method
        
    end
    
    
end

