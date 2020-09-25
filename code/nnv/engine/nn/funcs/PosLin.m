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
        function S = stepReach(varargin)
            % @I: single star set input
            % @index: index of the neuron performing stepPosLin
            % @xmin: minimum of x[index]
            % @xmax: maximum of x[index]
            % @S: star output set
            
            % author: Dung Tran
            % date: 27/2/2019
            % update: 7/16/2020: add lp_solver option
            
            switch nargin
                case 2 
                    I = varargin{1};
                    index = varargin{2};
                    lp_solver = 'linprog';
                case 3
                    I = varargin{1};
                    index = varargin{2};
                    lp_solver = varargin{3};
                otherwise
                    error('Invalid number of input arguments, should be 2 or 3');
            end
            
            if ~isa(I, 'Star')
                error('Input is not a star set');
            end
            
            xmin = I.getMin(index, lp_solver);
                       
            if xmin >= 0
                S = I; 
            else
                xmax = I.getMax(index, lp_solver);
                
                if xmax <= 0
                
                    V1 = I.V;
                    V1(index, :) = 0;
                    if ~isempty(I.Z)
                        c = I.Z.c;
                        c(index) = 0;
                        V = I.Z.V;
                        V(index, :) = 0;
                        new_Z = Zono(c, V); % update outer-zono
                    else
                        new_Z = [];
                    end
                    S = Star(V1, I.C, I.d, I.predicate_lb, I.predicate_ub, new_Z);

                else
                    % S1 = I && x[index] < 0 
                    c = I.V(index, 1);
                    V = I.V(index, 2:I.nVar + 1); 
                    new_C = vertcat(I.C, V);
                    new_d = vertcat(I.d, -c);                
                    new_V = I.V;
                    new_V(index, :) = zeros(1, I.nVar + 1);

                    % update outer-zono
                    if ~isempty(I.Z)
                        c1 = I.Z.c;
                        c1(index) = 0;
                        V1 = I.Z.V;
                        V1(index, :) = 0;
                        new_Z = Zono(c1, V1);
                    else
                        new_Z = [];
                    end
                    S1 = Star(new_V, new_C, new_d, I.predicate_lb, I.predicate_ub, new_Z);

                    % S2 = I && x[index] >= 0
                    new_C = vertcat(I.C, -V);
                    new_d = vertcat(I.d, c);

                    S2 = Star(I.V, new_C, new_d, I.predicate_lb, I.predicate_ub, I.Z);

                   S = [S1 S2];

                end
            end

        end
        
        % new stepReach method, compute reachable set for a single step
        % minimize the number of LP using simulation 
        function S = stepReach2(I, index)
            % @I: single star set input
            % @index: index of the neuron performing stepPosLin
            % @xmin: minimum of x[index]
            % @xmax: maximum of x[index]
            % @S: star output set
            
            % author: Dung Tran
            % date: 27/2/2019
            
            if ~isa(I, 'Star')
                error('Input is not a star set');
            end
            
            x1 = I.V(index, 1) + I.V(index, 2:I.nVar + 1)*I.predicate_lb;
            x2 = I.V(index, 1) + I.V(index, 2:I.nVar + 1)*I.predicate_ub;
            
            if x1*x2 < 0
                % S1 = I && x[index] < 0 
                    c = I.V(index, 1);
                    V = I.V(index, 2:I.nVar + 1); 
                    new_C = vertcat(I.C, V);
                    new_d = vertcat(I.d, -c);                
                    new_V = I.V;
                    new_V(index, :) = zeros(1, I.nVar + 1);

                    % update outer-zono
                    if ~isempty(I.Z)
                        c1 = I.Z.c;
                        c1(index) = 0;
                        V1 = I.Z.V;
                        V1(index, :) = 0;
                        new_Z = Zono(c1, V1);
                    else
                        new_Z = [];
                    end
                    S1 = Star(new_V, new_C, new_d, I.predicate_lb, I.predicate_ub, new_Z);

                    % S2 = I && x[index] >= 0
                    new_C = vertcat(I.C, -V);
                    new_d = vertcat(I.d, c);
                    S2 = Star(I.V, new_C, new_d, I.predicate_lb, I.predicate_ub, I.Z);
                    
                    S = [S1 S2];

            else
                
                if (x1 < 0 && x2 < 0)
                    xmax = I.getMax(index);
                    if xmax <= 0

                        V1 = I.V;
                        V1(index, :) = 0;
                        if ~isempty(I.Z)
                            c = I.Z.c;
                            c(index) = 0;
                            V = I.Z.V;
                            V(index, :) = 0;
                            new_Z = Zono(c, V); % update outer-zono
                        else
                            new_Z = [];
                        end
                        S = Star(V1, I.C, I.d, I.predicate_lb, I.predicate_ub, new_Z);

                    else
                        % S1 = I && x[index] < 0 
                        c = I.V(index, 1);
                        V = I.V(index, 2:I.nVar + 1); 
                        new_C = vertcat(I.C, V);
                        new_d = vertcat(I.d, -c);                
                        new_V = I.V;
                        new_V(index, :) = zeros(1, I.nVar + 1);

                        % update outer-zono
                        if ~isempty(I.Z)
                            c1 = I.Z.c;
                            c1(index) = 0;
                            V1 = I.Z.V;
                            V1(index, :) = 0;
                            new_Z = Zono(c1, V1);
                        else
                            new_Z = [];
                        end
                        S1 = Star(new_V, new_C, new_d, I.predicate_lb, I.predicate_ub, new_Z);

                        % S2 = I && x[index] >= 0
                        new_C = vertcat(I.C, -V);
                        new_d = vertcat(I.d, c);

                        S2 = Star(I.V, new_C, new_d, I.predicate_lb, I.predicate_ub, I.Z);

                       S = [S1 S2];

                    end
                else
                    
                    xmin = I.getMin(index);
                       
                    if xmin >= 0
                        S = I; 
                    else
                        
                        % S1 = I && x[index] < 0 
                        c = I.V(index, 1);
                        V = I.V(index, 2:I.nVar + 1); 
                        new_C = vertcat(I.C, V);
                        new_d = vertcat(I.d, -c);                
                        new_V = I.V;
                        new_V(index, :) = zeros(1, I.nVar + 1);

                        % update outer-zono
                        if ~isempty(I.Z)
                            c1 = I.Z.c;
                            c1(index) = 0;
                            V1 = I.Z.V;
                            V1(index, :) = 0;
                            new_Z = Zono(c1, V1);
                        else
                            new_Z = [];
                        end
                        S1 = Star(new_V, new_C, new_d, I.predicate_lb, I.predicate_ub, new_Z);

                        % S2 = I && x[index] >= 0
                        new_C = vertcat(I.C, -V);
                        new_d = vertcat(I.d, c);

                        S2 = Star(I.V, new_C, new_d, I.predicate_lb, I.predicate_ub, I.Z);

                       S = [S1 S2];

                    end
                                       
                end
                
            end

        end
        
        
        % stepReach with multiple inputs
        function S = stepReachMultipleInputs(varargin)
            % @I: an array of stars
            % @index: index where stepReach is performed
            % @option: = 'parallel' use parallel computing
            %          = not declare -> don't use parallel computing
            
            % author: Dung Tran
            % date: 27/2/2019
            % update: 4/2/2020
            %       : 7/16/2020: add lp_solver option
            
            switch nargin
                case 3
                    I = varargin{1};
                    index = varargin{2};
                    option = varargin{3};
                    lp_solver = 'linprog';
                case 4
                    I = varargin{1};
                    index = varargin{2};
                    option = varargin{3};
                    lp_solver = varargin{4};
                otherwise
                    error('Invalid number of input arguments');
            end
            
            p = length(I);
            S = [];
            
            if isempty(option)
                
                for i=1:p
                    S =[S, PosLin.stepReach(I(i), index, lp_solver)];
                end
                
            elseif strcmp(option, 'parallel')
                
                parfor i=1:p
                    S =[S, PosLin.stepReach(I(i), index, lp_solver)];
                end
                
            else
                error('Unknown option');
            end
            
            
        end      
        
                
        
        % exact reachability analysis using star
        function S = reach_star_exact(varargin)
            % @I: star input sets
            % @option: = 'parallel' using parallel computing
            %          = '[]'    do not use parallel computing
            
            % author: Dung Tran
            % date: 3/16/2019
            % update: 7/15/2020: add display option
            %         7/16/2020: add lp_solver option
            
            switch nargin
                case 2
                    I = varargin{1};
                    option = varargin{2};
                    dis_opt = [];
                    lp_solver = 'linprog';
                case 3
                    I = varargin{1};
                    option = varargin{2};
                    dis_opt = varargin{3};
                    lp_solver = 'linprog';
                case 4
                    I = varargin{1};
                    option = varargin{2};
                    dis_opt = varargin{3};
                    lp_solver = varargin{4};
                otherwise
                    error('Invalid number of input arguments');
            end
            
             if ~isempty(I)
                            
                [lb, ub] = I.estimateRanges;
                
                if isempty(lb) || isempty(ub)
                    S = [];
                else
                    map = find(ub <= 0); % computation map
                    V = I.V;
                    V(map, :) = 0;
                    % update outer-zono
                    if ~isempty(I.Z)
                        c1 = I.Z.c;
                        c1(map, :) = 0;
                        V1 = I.Z.V;
                        V1(map, :) = 0;
                        new_Z = Zono(c1, V1);
                    else
                        new_Z = [];
                    end
                    
                    In = Star(V, I.C, I.d, I.predicate_lb, I.predicate_ub, new_Z);                    
                    map = find(lb < 0 & ub > 0);
                    m = length(map);                    
                    for i=1:m
                        if strcmp(dis_opt, 'display')
                            fprintf('\nPerforming exact PosLin_%d operation using Star', map(i));
                        end
                        In = PosLin.stepReachMultipleInputs(In, map(i), option, lp_solver);
                    end               
                    S = In;
                end
                
            else
                S = [];
            end
            
        end
        
        
%         % a new implementation of exact-star method
%         % generate multiple Stars at the same time
%         function S = reach_star_exact2(I, ~)
%             % @I: star input sets
%             % @option: = 'parallel' using parallel computing
%             %          = '[]'    do not use parallel computing
%             
%             % author: Dung Tran
%             % date: 6/16/2020
%             
%             
%              if ~isempty(I)
%                             
%                 [lb, ub] = I.estimateRanges; 
%                 map1 = find(ub <= 0); % reset map
%                 map2 = find(lb < 0 & ub >0); % coarse split map
% 
%                 x1 = I.V(:, 1) + I.V(:,2:I.nVar+1) * I.predicate_lb; % first test point
%                 x2 = I.V(:, 1) + I.V(:,2:I.nVar+1) * I.predicate_ub; % second test point
%                 map3 = find(x1 < 0 & x2 > 0);
%                 map4 = find(x1 > 0 & x2 < 0);
% 
%                 map5 = [map3; map4]; % obvious split map
%                 map6 = setdiff(map2, map5); % a reduced map to find optimized ranges
% 
%                 xmax = I.getMaxs(map6,'no-show-progess');
%                 map7 = find(xmax <= 0); 
%                 map8 = map6(map7(:)); % the index that has xmax <= 0
%                 map9 = [map1; map8]; % the final reset map that contains all neurons whose values <= 0
% 
%                 map71 = setdiff(map6, map8); % the index that has xmax > 0
%                 xmin = I.getMins(map71, 'no-show-progess'); % minimum of the index that has xmax <= 0
% 
%                 map10 = find(xmin < 0);
%                 map11 = map71(map10(:)); % the neuron indexes that has xmin <0 & xmax > 0
% 
%                 map12 = [map5; map11]; % the final slit map containing all neurons that can split, i.e. values contain 0
% 
%                 %n_LP = length(map6) + length(map8);
%                 %fprintf('\nNeglected LP: %d/%d = %.3f5%', I.dim - n_LP, I.dim, ((I.dim-n_LP)*100)/I.dim);
% 
%                 V1 = I.V; 
%                 V1(map9,:) = 0; % reset to zero all neurons whose values <= 0              
% 
%                 N = length(map12); % there are N neurons to split
%                 M = 0:2^N-1; % there are 2^N new star sets
%                 A = de2bi(M)'; % Nx2^N matrix to memorize where the plit is, 0 <-> x <0, 1 <-> x > 0
% 
%                 % construct 2^N star sets at the same time
%                 G = length(M);
%                 S = [];
%                 for i=1:G
%                     V2 = V1;
%                     V2(map12, :) = A(:,i).*V1(map12, :);
%                     C21 = -A(:,i).*V1(map12, 2:I.nVar+1);% x >0
%                     d21 = A(:,i).*V1(map12, 1);
%                     C22 = (1-A(:,i)).*V1(map12, 2:I.nVar+1); % x < 0
%                     d22 = -(1-A(:,i)).*V1(map12,1);
%                     C2 = C21 + C22;
%                     d2 = d21 + d22;
%                     I_new = Star(V2, [I.C; C2], [I.d; d2], I.predicate_lb, I.predicate_ub);
%                     if ~I_new.isEmptySet
%                         S = [S I_new];
%                     end
%                 end
% 
%                
%             else
%                 S = [];
%             end
%             
%         end
        
        
        % exact reachability analysis using star
        function S = reach_star_exact_multipleInputs(varargin)
            % @In: star input sets
            % @option: = 'parallel' using parallel computing
            %          = '[]'    do not use parallel computing
            
            % author: Dung Tran
            % date: 7/25/2019
            % update: 7/16/2020: add display option + lp_solver option
            
            switch nargin
                case 2
                    In = varargin{1};
                    option = varargin{2};
                    dis_opt = [];
                    lp_solver = 'linprog';
                case 3
                    In = varargin{1};
                    option = varargin{2};
                    dis_opt = varargin{3};
                    lp_solver = 'linprog';
                case 4
                    In = varargin{1};
                    option = varargin{2};
                    dis_opt = varargin{3};
                    lp_solver = varargin{4};
                otherwise
                    error('Invalid number of input arguments, should be 2, 3 or 4');
            end
            
            
             n = length(In);
             S = [];
             if strcmp(option, 'parallel')
                 parfor i=1:n
                     S = [S PosLin.reach_star_exact(In(i), [], dis_opt, lp_solver)];
                 end
             elseif isempty(option) || strcmp(option, 'single')
                 for i=1:n
                     S = [S PosLin.reach_star_exact(In(i), [], dis_opt, lp_solver)];
                 end
             else
                 error('Unknown computation option');
             end
        end
              
        
        
        % step reach approximation using star
        function S = stepReachStarApprox(I, index)
            % @I: star set input
            % @index: index of the neuron performing stepReach
            % @S: star output

            % author: Dung Tran
            % date: 7/22/2019
                       
            
            if ~isa(I, 'Star')
                error('Input is not a star');
            end
                        
            lb = I.getMin(index);
              
            if lb > 0
                S = I;
            else
                ub = I.getMax(index);
                if ub <= 0
                    V = I.V;
                    V(index, :) = 0;
                    if ~isempty(I.Z)
                        c1= I.Z.c;
                        c1(index) = 0;
                        V1 = I.Z.V;
                        V1(index, :) = 0;
                        new_Z = Zono(c1, V1); % update outer-zono
                    else
                        new_Z = [];
                    end
                    S = Star(V, I.C, I.d, I.predicate_lb, I.predicate_ub, new_Z);
                else
                    fprintf('\nAdd a new predicate variables at index = %d', index);
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
                    new_predicate_lb = [I.predicate_lb; 0];                
                    new_predicate_ub = [I.predicate_ub; ub];

                    % update outer-zono

                    lamda = ub/(ub -lb);
                    mu = -0.5*ub*lb/(ub - lb);
                    if ~isempty(I.Z)
                        c = I.Z.c; 
                        c(index) = lamda * c(index) + mu;
                        V = I.Z.V;
                        V(index, :) = lamda * V(index, :);
                        I1 = zeros(I.dim,1);
                        I1(index) = mu;
                        V = [V I1];              
                        new_Z = Zono(c, V);
                    else
                        new_Z = [];
                    end

                    S = Star(new_V, new_C, new_d, new_predicate_lb, new_predicate_ub, new_Z);
                end
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

            if isempty(I)
                S = [];
            else
                [lb, ub] = I.estimateRanges;
                if isempty(lb) || isempty(ub)
                    S = [];
                else
                    map = find(ub <= 0); % computation map
                    V = I.V;
                    V(map, :) = 0;
                    if ~isempty(I.Z)
                        c1 = I.Z.c;
                        c1(map, :) = 0;
                        V1 = I.Z.V;
                        V1(map, :) = 0;
                        new_Z = Zono(c1, V1);
                    else
                        new_Z = [];
                    end
                    In = Star(V, I.C, I.d, I.predicate_lb, I.predicate_ub, new_Z);                    
                    map = find(lb < 0 & ub > 0);
                    m = length(map); 
                    for i=1:m
                        fprintf('\nPerforming approximate PosLin_%d operation using Star', map(i));
                        In = PosLin.stepReachStarApprox(In, map(i));
                    end
                    S = In;
             
                end
            end

        end
        
        
        % step reach approximation using star
        function S = multipleStepReachStarApprox_at_one(I, index, lb, ub)
            % @I: star set input
            % @index: index of the neurons performing stepReach
            % @lb:lower bound of x[index]
            % @ub: upper bound of x[index]

            % author: Dung Tran
            % date: 6/11/2020
                       
            if ~isa(I, 'Star')
                error('Input is not a star');
            end
            
            N = I.dim; 
            m = length(index); % number of neurons involved (number of new predicate variables introduced)
            
            % construct new basis array
            V1 = I.V; % originial basis array
            V1(index, :) = 0;
            V2 = zeros(N, m); % basis array for new predicates
            for i=1:m
                V2(index(i), i) = 1;
            end
            new_V = [V1 V2]; % new basis for over-approximate star set
            
            % construct new constraints on new predicate variables
            
            % case 0: keep the old constraints on the old predicate
            % variables
            
            n = I.nVar; % number of old predicate variables
            C0 = [I.C zeros(size(I.C, 1), m)];
            d0 = I.d; 
            
            %case 1: y(index) >= 0           
            C1 = [zeros(m, n) -eye(m)];
            d1 = zeros(m, 1);
            
            %case 2: y(index) >= x(index)
            C2 = [I.V(index, 2:n+1) -V2(index, 1:m)];
            d2 = -I.V(index, 1);
            
            %case 3: y(index) <= (ub/(ub - lb))*(x-lb)
            a = ub./(ub - lb); % divide element-wise
            b = a.*lb; % multiply element-wise
            C3 = [-a.*I.V(index, 2:n+1) V2(index, 1:m)];
            d3 = a.*I.V(index, 1)-b;
            
            new_C = [C0; C1; C2; C3];
            new_d = [d0; d1; d2; d3];
            
            new_pred_lb = [I.predicate_lb; zeros(m,1)];
            new_pred_ub = [I.predicate_ub; ub];
            
            S = Star(new_V, new_C, new_d, new_pred_lb, new_pred_ub);
            
        end
        
        
        % more efficient method by doing multiple stepReach at one time
        % over-approximate reachability analysis using Star
        function S = reach_star_approx2(varargin)
            % @I: star input set
            % @option: 'parallel' or single
            % @S: star output set

            % author: Dung Tran
            % date: 4/3/2019
            % update: 7/13/2020 : getMax parallel
            % update: 7/15/2020 : add display option
            %         7/16/2020: add lp_solver option
            
            switch nargin
                case 1
                    I = varargin{1};
                    option = 'single';
                    dis_opt = [];
                    lp_solver = 'linprog';
                case 2
                    I = varargin{1};
                    option = varargin{2};
                    dis_opt = [];
                    lp_solver = 'linprog';
                case 3
                    I = varargin{1};
                    option = varargin{2};
                    dis_opt = varargin{3};
                    lp_solver = 'linprog';
                case 4
                    I = varargin{1};
                    option = varargin{2};
                    dis_opt = varargin{3};
                    lp_solver = varargin{4};
                otherwise
                    error('Invalid number of input arguments, should be 1, 2, 3, or 4');
            end

            if ~isa(I, 'Star')
                error('Input is not a star');
            end

            if isempty(I)
                S = [];
            else
                [lb, ub] = I.estimateRanges;
                if isempty(lb) || isempty(ub)
                    S = [];
                else
                    % find all indexes having ub <= 0, then reset the
                    % values of the elements corresponding to these indexes to 0
                    if strcmp(dis_opt, 'display')
                        fprintf('\nFinding all neurons (in %d neurons) with ub <= 0...', length(ub));
                    end
                    map1 = find(ub <= 0); % computation map
                    if strcmp(dis_opt, 'display')
                        fprintf('\n%d neurons with ub <= 0 are found by estimating ranges', length(map1));
                    end

                    map2 = find(lb < 0 & ub > 0);
                    if strcmp(dis_opt, 'display')
                        fprintf('\nFinding neurons (in %d neurons) with ub <= 0 by optimizing ranges: ', length(map2));
                    end
                    xmax = I.getMaxs(map2, option, dis_opt, lp_solver);
                    map3 = find(xmax <= 0);
                    if strcmp(dis_opt, 'display')
                        fprintf('\n%d neurons (in %d neurons) with ub <= 0 are found by optimizing ranges', length(map3), length(map2));
                    end
                    n = length(map3);
                    map4 = zeros(n,1);
                    for i=1:n
                        map4(i) = map2(map3(i));
                    end
                    map11 = [map1; map4];
                    In = I.resetRow(map11); % reset to zero at the element having ub <= 0
                    if strcmp(dis_opt, 'display')
                        fprintf('\n(%d+%d =%d)/%d neurons have ub <= 0', length(map1), length(map3), length(map11), length(ub));
                    end

                    % find all indexes that have lb < 0 & ub > 0, then
                    % apply the over-approximation rule for ReLU
                    if strcmp(dis_opt, 'display')
                        fprintf("\nFinding all neurons (in %d neurons) with lb < 0 & ub >0: ", length(ub));
                    end
                    map5 = find(xmax > 0);
                    map6 = map2(map5(:)); % all indexes having ub > 0
                    xmax1 = xmax(map5(:)); % upper bound of all neurons having ub > 0

                    xmin = I.getMins(map6, option, dis_opt, lp_solver); 
                    map7 = find(xmin < 0); 
                    map8 = map6(map7(:)); % all indexes having lb < 0 & ub > 0
                    lb1 = xmin(map7(:));  % lower bound of all indexes having lb < 0 & ub > 0
                    ub1 = xmax1(map7(:)); % upper bound of all neurons having lb < 0 & ub > 0
                    
                    if strcmp(dis_opt, 'display')
                        fprintf('\n%d/%d neurons have lb < 0 & ub > 0', length(map8), length(ub));
                        fprintf('\nConstruct new star set, %d new predicate variables are introduced', length(map8));
                    end
                    S = PosLin.multipleStepReachStarApprox_at_one(In, map8, lb1, ub1); % one-shot approximation
                end
            end

        end        
        
        % a relaxed star-approx method
        function S = reach_relaxed_star_approx(varargin)
            % @I: star input set
            % @relaxFactor: a relaxFactor
            % @S: star output set

            % author: Dung Tran
            % date: 6/26/2020
            % update: 7/15/2020 add display option
            %       : 7/16/2020 add lp_solver option
            
            switch nargin
                case 2
                    I = varargin{1};
                    relaxFactor = varargin{2};
                    option = 'single';
                    dis_opt = [];
                    lp_solver = 'linprog';
                case 3
                    I = varargin{1};
                    relaxFactor = varargin{2};
                    option = varargin{3};
                    dis_opt = [];
                    lp_solver = 'linprog';
                case 4
                    I = varargin{1};
                    relaxFactor = varargin{2};
                    option = varargin{3};
                    dis_opt = varargin{4};
                    lp_solver = 'linprog';
                case 5
                    I = varargin{1};
                    relaxFactor = varargin{2};
                    option = varargin{3};
                    dis_opt = varargin{4};
                    lp_solver = varargin{5};
                otherwise
                    error('Invalid number of input arguments, should be 2, 3, 4 or 5');
            end

            if ~isa(I, 'Star')
                error('Input is not a star');
            end
            if relaxFactor < 0 || relaxFactor > 1
                error('Invalid relax factor');
            end

            if isempty(I)
                S = [];
            else
                [lb, ub] = I.estimateRanges;
                if isempty(lb) || isempty(ub)
                    S = [];
                else
                    
                    % find all indexes having ub <= 0, then reset the
                    % values of the elements corresponding to these indexes to 0
                    if strcmp(dis_opt, 'display')
                        fprintf('\nFinding all neurons (in %d neurons) with ub <= 0...', length(ub));
                    end
                    map1 = find(ub <= 0); % computation map
                    if strcmp(dis_opt, 'display')
                        fprintf('\n%d neurons with ub <= 0 are found by estimating ranges', length(map1));
                    end
                    map2 = find(lb < 0 & ub > 0);
                    n1  = round((1-relaxFactor)*length(map2)); % number of LP need to solve
                    if strcmp(dis_opt, 'display')
                        fprintf('\nFinding neurons (in (1-%.3f) x %d neurons = %d) with ub <= 0 by optimizing ranges, i.e. relaxing %2.2f%%: ', relaxFactor, length(map2), n1, 100*relaxFactor);                 
                    end
                    [~,midx] = sort(ub(map2)-lb(map2), 'descend');
                    map21 = map2(midx(1:n1)); % neurons with optimized ranged
                    map22 = map2(midx(n1+1:length(map2))); % neurons without optimized ranges
                    lb1 = lb(map22);
                    ub1 = ub(map22); 
                    
                    xmax = I.getMaxs(map21, option, dis_opt, lp_solver); 
                    map3 = find(xmax <= 0);
                    if strcmp(dis_opt, 'display')
                        fprintf('\n%d neurons (in %d neurons) with ub <= 0 are found by optimizing ranges', length(map3), length(map21));
                    end
                    n = length(map3);
                    map4 = zeros(n,1);
                    for i=1:n
                        map4(i) = map2(map3(i));
                    end
                    map11 = [map1; map4];
                    In = I.resetRow(map11); % reset to zero at the element having ub <= 0
                    if strcmp(dis_opt, 'display')
                        fprintf('\n(%d+%d =%d)/%d neurons have ub <= 0', length(map1), length(map3), length(map11), length(ub));
                    end

                    % find all indexes that have lb < 0 & ub > 0, then
                    % apply the over-approximation rule for ReLU
                    if strcmp(dis_opt, 'display')
                        fprintf("\nFinding neurons (in %d neurons) with lb < 0 & ub >0: ", length(map21));
                    end
                    map5 = find(xmax > 0);
                    map6 = map21(map5(:)); % all indexes having ub > 0
                    xmax1 = xmax(map5(:)); % upper bound of all neurons having ub > 0

                    xmin = I.getMins(map6, option, dis_opt, lp_solver); 
                    map7 = find(xmin < 0); 
                    map8 = map6(map7(:)); % all indexes having lb < 0 & ub > 0
                    lb2 = xmin(map7(:));  % lower bound of all indexes having lb < 0 & ub > 0
                    ub2 = xmax1(map7(:)); % upper bound of all neurons having lb < 0 & ub > 0
                    
                    map9 = [map22; map8];
                    lb3 = [lb1; lb2];
                    ub3 = [ub1; ub2];
                    if strcmp(dis_opt, 'display')
                        fprintf('\n%d/%d neurons have lb < 0 & ub > 0', length(map9), length(ub));
                        fprintf('\nConstruct new star set, %d new predicate variables are introduced', length(map9));
                    end
                    S = PosLin.multipleStepReachStarApprox_at_one(In, map9, lb3, ub3); % one-shot approximation
                end
            end

        end


    end
    
    
    methods(Static) % reachability analysis using Polyhedron
        
        % step reach set y = ReLU(x)
        function R = stepReach_Polyhedron(I, index, xmin, xmax)
            % @I : input set, a polyhedron
            % @i : index of current x[index] of current step
            % @xmin: min of x[index]
            % @xmax: max of x[index]
           
            I.normalize;
            dim = I.Dim;
            if xmin >= 0
                R = I; 
            elseif xmax < 0 
                Im = eye(dim);
                Im(index, index) = 0;
                R = I.affineMap(Im, 'vrep');
            elseif xmin < 0 && xmax >= 0
                
                Z1 = zeros(1, dim);
                Z1(1, index) = 1;
                Z2 = zeros(1, dim);
                Z2(1, index) = -1;

                A1 = vertcat(I.A, Z1);
                A2 = vertcat(I.A, Z2);
                b  = vertcat(I.b, [0]);
                R1 = Polyhedron('A', A1, 'b', b, 'Ae', I.Ae, 'be', I.be);
                R2 = Polyhedron('A', A2, 'b', b, 'Ae', I.Ae, 'be', I.be);
                
                Im = eye(dim);
                Im(index, index) = 0;
                R1 = R1.affineMap(Im, 'vrep');
                if R1.isEmptySet 
                    if R2.isEmptySet
                        R = [];
                    else
                        
                        R = R2;
                    end
                else
                    if R2.isEmptySet
                        R = R1;
                    else
                        R = [R1 R2];
                    end
                    
                end               
                
            end
            
        end
        
        
        % stepReach for multiple Input Sets 
        function R = stepReachMultipleInputs_Polyhedron(varargin)
            % @I: an array of input sets which are polyhedra
            % @index: index of current x[index] of current step
            % @xmin: min value of x[index]
            % @xmax: max value of x[index]
            % @option: = 'exact' -> compute an exact reach set
            %          = 'approx' -> compute an over-approximate reach set
            
            switch nargin
                
                case 5
                    I = varargin{1};
                    index = varargin{2};
                    xmin = varargin{3};
                    xmax = varargin{4};
                    option = varargin{5};
                
                otherwise
                    error('Invalid number of input arguments (should be 5)');
            end
            
            
            
            p = length(I);
            R = [];
            
            if isempty(option)
                
                for i=1:p
                    R =[R, PosLin.stepReach_Polyhedron(I(i), index, xmin, xmax)];
                end
                
            elseif strcmp(option, 'parallel')
                
                parfor i=1:p
                    R =[R, PosLin.stepReach_Polyhedron(I(i), index, xmin, xmax)];
                end
                
            else
                error('Unknown option');
            end
            
            
                     
        end
        
        
        % exact reachability analysis using Polyhedron
        function R = reach_polyhedron_exact(varargin)
            
                        % @I: star input sets
            % @option: = 'parallel' using parallel computing
            %          = ''    do not use parallel computing
            
            % author: Dung Tran
            % date: 3/16/2019
            % update: 7/15/2020 : add display option
            
            switch nargin
                case 2
                    I = varargin{1};
                    option = varargin{2};
                    dis_opt = [];
                case 3
                    I = varargin{1};
                    option = varargin{2};
                    dis_opt = varargin{3};
                otherwise
                    error('Invalid number of input arguments, should be 2 or 3');
            end
            
             if ~isempty(I)
                if isa(I, 'Polyhedron')            
                    I.outerApprox; % find bounds of I state vector
                    lb = I.Internal.lb; % min-vec of x vector
                    ub = I.Internal.ub; % max-vec of x vector
                else
                    error('Input set is not a Polyhedron');
                end
                
                if isempty(lb) || isempty(ub)
                    R = [];
                else
                    map = find(lb < 0); % computation map
                    m = size(map, 1); % number of stepReach operations needs to be executed
                    In = I;
                    for i=1:m
                        if strcmp(dis_opt, 'display')
                            fprintf('\nPerforming exact PosLin_%d operation using Polyhedron', map(i));
                        end
                        In = PosLin.stepReachMultipleInputs_Polyhedron(In, map(i), lb(map(i)), ub(map(i)), option);
                    end               
                    R = In;
                end
                
            else
                R = [];
            end
                 
        end
             
    end


    methods(Static) % over-approximate reachability analysis use zonotope
        
        % step over-approximate reachability analysis using zonotope
        function Z = stepReachZonoApprox(I, index, lb, ub)
            % @I: zonotope input set
            % @lb: lower bound of input at specific neuron i
            % @ub: lower bound of input at specfic neuron i
            % @index: index of the neuron we want to perform stepReach
            % @Z: zonotope output set
            
            % author: Dung Tran
            % date: 5/3/2019
            % update: 1/4/2020:
            %     update reason: change getBounds method for Zonotope
        
            % reference: Fast and Effective Robustness Ceritification,
            % Gagandeep Singh, NIPS 2018
           
            if ~isa(I, 'Zono')
                error('Input is not a Zonotope');
            end
            
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
        function Z = reach_zono_approx(varargin)
            % @I: zonotope input
            % @Z: zonotope output
            
            % author: Dung Tran
            % date: 5/3/2019
            
            % reference: Fast and Effective Robustness Ceritification,
            % Gagandeep Singh, NIPS 2018
            
            % update: 7/15/2020 : add display option
            
            switch nargin
                case 1
                    I = varargin{1};
                    dis_opt = [];
                case 2
                    I = varargin{1};
                    dis_opt = varargin{2};
                otherwise
                    error('Invalid number of input arguments, should be 1 or 2');
            end
            
            if ~isa(I, 'Zono')
                error('Input is not a Zonotope');
            end
                      
            In = I;
            [lb, ub] = I.getBounds;
            for i=1:I.dim
                if strcmp(dis_opt, 'display')
                    fprintf('\nPerforming approximate PosLin_%d operation using Zonotope', i);
                end
                In = PosLin.stepReachZonoApprox(In, i, lb(i), ub(i));
            end
            Z = In;
                       
        end
        
        
    end
    
    
    methods(Static) % reachability analysis using abstract-domain
        
        % step over-approximate reachability analysis using abstract-domain
        % we use star set to represent abstract-domain
        function S = stepReachAbstractDomain(varargin)
            % @I: star-input set
            % @index: index of neuron performing stepReach
            % @S: star output set represent abstract-domain of the output
            % set
            
            % author: Dung Tran
            % date: 5/3/2019
        
            % reference: An Abstract Domain for Certifying Neural Networks,
            % Gagandeep Singh, POPL 2019
            
            
            switch nargin
                
                case 4
                    
                    I = varargin{1};
                    index = varargin{2};
                    lb = varargin{3};
                    ub = varargin{4};
                
                case 2
                    I = varargin{1};
                    index = varargin{2};
                    %[lb, ub] = I.getRange(index); % our improved approach
                    [lb, ub] = I.estimateRange(index); % originial DeepPoly approach use estimated range
                otherwise
                    error('Invalid number of input arguments (should be 2 or 4)');
            end
            
            if ~isa(I, 'Star')
                error('Input is not a Star');
            end
                          
            if lb > 0
                S = Star(I.V, I.C, I.d, I.predicate_lb, I.predicate_ub);
            elseif ub <= 0
                V = I.V;
                V(index, :) = zeros(1, I.nVar + 1);
                S = Star(V, I.C, I.d, I.predicate_lb, I.predicate_ub);
            elseif lb <= 0 && ub > 0
                
                S1 = ub*(ub-lb)/2; % area of the first candidate abstract-domain
                S2 = -lb*(ub-lb)/2; % area of the second candidate abstract-domain  
                
                n = I.nVar + 1;
                                
                % constraint 1: y[index] = ReLU(x[index]) >= 0
                C1 = zeros(1, n);
                C1(n) = -1; 
                d1 = 0;
                
                
                % constraint 2: y[index] = ReLU(x[index]) >= x[index]
                C2 = [I.V(index, 2:n) -1];
                d2 = -I.V(index, 1);
                    
                % constraint 3: y[index] <= ub(x[index] -lb)/(ub - lb)
                C3 = [-(ub/(ub-lb))*I.V(index, 2:n) 1];
                d3 = -(ub*lb/(ub-lb)) + ub*I.V(index, 1)/(ub-lb);
                               
                m = size(I.C, 1);
                C0 = [I.C zeros(m, 1)];
                d0 = I.d;
                new_V = [I.V zeros(I.dim, 1)];
                new_V(index, :) = zeros(1, n+1);
                new_V(index, n+1) = 1;
                
                if S1 < S2
                    % get first cadidate as resulted abstract-domain
                    new_C = [C0; C1; C3];
                    new_d = [d0; d1; d3];
                    new_pred_lb = [I.predicate_lb; 0];
                    new_pred_ub = [I.predicate_ub; ub];
                    
                else
                    % choose the second candidate as the abstract-domain                                      
                    new_C = [C0; C2; C3];
                    new_d = [d0; d2; d3];
                    new_pred_lb = [I.predicate_lb; lb];
                    new_pred_ub = [I.predicate_ub; ub];
                                        
                end
		display(new_C)
                S = Star(new_V, new_C, new_d, new_pred_lb, new_pred_ub);
                
            end
                       
        end
        
        
        % over-approximate reachability analysis using abstract-domain
        function S = reach_abstract_domain(varargin)
            % @I: star input set
            % @S: star output set

            % author: Dung Tran
            % date: 4/3/2019

            % update: 7/15/2020 : add display option
            
            switch nargin
                case 1
                    I = varargin{1};
                    dis_opt = [];
                case 2
                    I = varargin{1};
                    dis_opt = varargin{2};
                otherwise
                    error('Invalid number of input arguments, should be 1 or 2');
            end
            
            if ~isa(I, 'Star')
                error('Input is not a star');
            end

            if isempty(I)
                S = [];
            else    
                [lb, ub] = I.estimateRanges;
                if isempty(lb) || isempty(ub)
                    S = [];
                else
                    map = find(ub <= 0); % computation map
                    V = I.V;
                    V(map, :) = 0;
                    In = Star(V, I.C, I.d, I.predicate_lb, I.predicate_ub);
                    map = find(lb <= 0 & ub > 0);
                    m = length(map); 
                    for i=1:m
                        if strcmp(dis_opt, 'display')
                            fprintf('\nPerforming approximate PosLin_%d operation using Abstract Domain', map(i));
                        end
                        In = PosLin.stepReachAbstractDomain(In, map(i), lb(map(i)), ub(map(i)));
                    end
                    S = In;
             
                end
            
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
            % update: 7/15/2020: add display option
            
            switch nargin
                
                case 6
                    I = varargin{1};
                    method = varargin{2};
                    option = varargin{3};
                    relaxFactor = varargin{4}; % used for aprox-star only
                    dis_opt = varargin{5}; % display option
                    lp_solver = varargin{6}; 
                
                case 5
                    I = varargin{1};
                    method = varargin{2};
                    option = varargin{3};
                    relaxFactor = varargin{4}; % used for aprox-star only
                    dis_opt = varargin{5}; % display option
                    lp_solver = 'linprog';
                
                case 4
                    I = varargin{1};
                    method = varargin{2};
                    option = varargin{3};
                    relaxFactor = varargin{4}; % used for aprox-star only
                    dis_opt = [];
                    lp_solver = 'linprog';
                                    
                case 3
                    I = varargin{1};
                    method = varargin{2};
                    option = varargin{3};
                    relaxFactor = 0; % used for aprox-star only
                    dis_opt = [];
                    lp_solver = 'linprog';
                case 2
                    I = varargin{1};
                    method = varargin{2};
                    option = 'parallel';
                    relaxFactor = 0; % used for aprox-star only
                    dis_opt = [];
                    lp_solver = 'linprog';
                case 1
                    I = varargin{1};
                    method = 'exact-star';
                    option = 'parallel';
                    relaxFactor = 0; % used for aprox-star only
                    dis_opt = [];
                    lp_solver = 'linprog';
                otherwise
                    error('Invalid number of input arguments (should be 1, 2, 3, 4, or 5)');
            end
            
            
            if strcmp(method, 'exact-star') % exact analysis using star
                
                R = PosLin.reach_star_exact_multipleInputs(I, option, dis_opt, lp_solver);
                
            elseif strcmp(method, 'exact-polyhedron') % exact analysis using polyhedron
                
                R = PosLin.reach_polyhedron_exact(I, option, dis_opt);
                
            elseif strcmp(method, 'approx-star')  % over-approximate analysis using star
                
                if relaxFactor == 0
                    R = PosLin.reach_star_approx2(I, option, dis_opt, lp_solver);
                else
                    R = PosLin.reach_relaxed_star_approx(I, relaxFactor, option, dis_opt, lp_solver);
                end
                
            elseif strcmp(method, 'approx-zono')  % over-approximate analysis using zonotope
                
                R = PosLin.reach_zono_approx(I, dis_opt);
                
            elseif strcmp(method, 'abs-dom')  % over-approximate analysis using abstract-domain
                
                R = PosLin.reach_abstract_domain(I, dis_opt);
                
            elseif strcmp(method, 'exact-face-latice') % exact analysis using face-latice
                fprintf('\nNNV have not yet support Exact Face-Latice Method');
            elseif strcmp(method, 'approx-face-latice') % over-approximate analysis using face-latice
                fprintf('\nNNV have not yet support Approximate Face-Latice Method');
            end
            
            
            
            
              
        end
        
        
    end
    
    
end

