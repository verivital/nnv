classdef LeakyReLU
    % LeakyReLU Class contains method for reachability analysis for Layer with leakyReLU activation function 
    
    % author: Dung Tran
    % date: 11/19/2020
    
    % TODO: why two reach_star_approx functions

    properties % no properties
        
    end
    
    methods(Static) % evaluate method and reachability analysis with stars 
        
        % evaluation
        function y = evaluate(x, gamma)
            y = x;
            y(find(y<0)) = gamma*y(find(y<0));
        end
        
        % stepReach method, compute reachable set for a single step
        function S = stepReach(varargin)
            % @I: single star set input
            % @index: index of the neuron performing stepLeakyReLU
            % @xmin: minimum of x[index]
            % @xmax: maximum of x[index]
            % @S: star output set
            
            % author: Dung Tran
            % date: 11/19/2020
            
            switch nargin
                case 3 
                    I = varargin{1};
                    index = varargin{2};
                    gamma = varargin{3};
                    lp_solver = 'linprog';
                case 4
                    I = varargin{1};
                    index = varargin{2};
                    gamma = varargin{3};
                    lp_solver = varargin{4};
                otherwise
                    error('Invalid number of input arguments, should be 3 or 4');
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
                    V1(index, :) = gamma*V1(index, :);
                    if ~isempty(I.Z)
                        c = I.Z.c;
                        c(index) = gamma*c(index);
                        V = I.Z.V;
                        V(index, :) = gamma*V(index,:);
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
                    new_V(index, :) = gamma*new_V(index, :);

                    % update outer-zono
                    if ~isempty(I.Z)
                        c1 = I.Z.c;
                        c1(index) = gamma*c1(index);
                        V1 = I.Z.V;
                        V1(index, :) = gamma*V1(index,:);
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
        
        % stepReach with multiple inputs
        function S = stepReachMultipleInputs(varargin)
            % @I: an array of stars
            % @index: index where stepReach is performed
            % @option: = 'parallel' use parallel computing
            %          = not declare -> don't use parallel computing
            
            % author: Dung Tran
            % date: 11/19/2020
            
            switch nargin
                case 3
                    I = varargin{1};
                    index = varargin{2};
                    gamma = varargin{3};
                    option = [];
                    lp_solver = 'linprog';
                case 4
                    I = varargin{1};
                    index = varargin{2};
                    gamma = varargin{3};
                    option = varargin{4};
                    lp_solver = 'linprog';
                case 5
                    I = varargin{1};
                    index = varargin{2};
                    gamma = varargin{3};
                    option = varargin{4};
                    lp_solver = varargin{5};
                otherwise
                    error('Invalid number of input arguments, should be 4 or 5');
            end
            
            p = length(I);
            S = [];
            
            if isempty(option)
                for i=1:p
                    S =[S LeakyReLU.stepReach(I(i), index, gamma, lp_solver)];
                end
            elseif strcmp(option, 'parallel')
                parfor i=1:p
                    S =[S LeakyReLU.stepReach(I(i), index, gamma, lp_solver)];
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
                case 3
                    I = varargin{1};
                    gamma = varargin{2};
                    option = varargin{3};
                    dis_opt = [];
                    lp_solver = 'linprog';
                case 4
                    I = varargin{1};
                    gamma = varargin{2};
                    option = varargin{3};
                    dis_opt = varargin{4};
                    lp_solver = 'linprog';
                case 5
                    I = varargin{1};
                    gamma = varargin{2};
                    option = varargin{3};
                    dis_opt = varargin{4};
                    lp_solver = varargin{5};
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
                    V(map, :) = gamma*V(map, :);
                    % update outer-zono
                    if ~isempty(I.Z)
                        c1 = I.Z.c;
                        c1(map, :) = gamma*c1(map, :);
                        V1 = I.Z.V;
                        V1(map, :) = gamma*V1(map, :);
                        new_Z = Zono(c1, V1);
                    else
                        new_Z = [];
                    end
                    
                    In = Star(V, I.C, I.d, I.predicate_lb, I.predicate_ub, new_Z);                    
                    map = find(lb < 0 & ub > 0);
                    m = length(map);                    
                    for i=1:m
                        if strcmp(dis_opt, 'display')
                            fprintf('\nPerforming exact LeakyReLU_%d operation using Star', map(i));
                        end
                        In = LeakyReLU.stepReachMultipleInputs(In, map(i), gamma, option, lp_solver);
                    end               
                    S = In;
                end
                
            else
                S = [];
            end
            
        end           
        
        % step reach approximation using star
        function S = stepReachStarApprox(varargin)
            % @I: star set input
            % @index: index of the neuron performing stepReach
            % @S: star output

            % author: Dung Tran
            % date: 7/22/2019
            % update: 11/23/2020
            
            switch nargin
                case 3
                    I = varargin{1};
                    index = varargin{2};
                    gamma = varargin{3};
                    lp_solver = 'linprog';
                case 4
                    I = varargin{1};
                    index = varargin{2};
                    gamma = varargin{3};
                    lp_solver = varargin{4};
                otherwise
                    error('Invalid number of input arguments, should be 3 or 4');
            end
            
            if ~isa(I, 'Star')
                error('Input is not a star');
            end
                        
            lb = I.getMin(index, lp_solver);
            if lb > 0
                S = I;
            else
                ub = I.getMax(index, lp_solver);
                if ub <= 0
                    V = I.V;
                    V(index, :) = gamma*V(index, :);
                    if ~isempty(I.Z)
                        c1= I.Z.c;
                        c1(index) = gamma*c1(index);
                        V1 = I.Z.V;
                        V1(index, :) = gamma*V1(index, :);
                        new_Z = Zono(c1, V1); % update outer-zono
                    else
                        new_Z = [];
                    end
                    S = Star(V, I.C, I.d, I.predicate_lb, I.predicate_ub, new_Z);
                else
                    n = I.nVar + 1;
                    % over-approximation constraints 
                    % constraint 1: y[index] = leakyReLU(x[index]) >= gamma*x[index]
                    V1 = I.V(index, 2:n);
                    c1 = I.V(index, 1);
                    C1 = [gamma*V1 -1];
                    d1 = -gamma*c1;
                    % constraint 2: y[index] >= x[index]
                    C2 = [V1 -1];
                    d2 = -c1;
                    % constraint 3: y[index] <= [(ub-gamma*lb)/(ub-lb)](x-lb) + gamma*lb
                    a = (ub-gamma*lb)/(ub-lb);
                    C3 = [-a*V1 1];
                    d3 = a*c1 - a*lb + gamma*lb;

                    m = size(I.C, 1);
                    C0 = [I.C zeros(m, 1)];
                    d0 = I.d;
                    new_C = [C0; C1; C2; C3];
                    new_d = [d0; d1; d2; d3];
                    new_V = [I.V zeros(I.dim, 1)];
                    new_V(index, :) = zeros(1, n+1);
                    new_V(index, n+1) = 1;              
                    new_predicate_lb = [I.predicate_lb; gamma*lb];                
                    new_predicate_ub = [I.predicate_ub; ub];

                    % update outer-zono
                    mu = (gamma*lb - a*lb)/2;
                    if ~isempty(I.Z)
                        c = I.Z.c; 
                        c(index) = a * c(index) + mu;
                        V = I.Z.V;
                        V(index, :) = a * V(index, :);
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
        function S = reach_star_approx(varargin)
            % @I: star input set
            % @S: star output set

            % author: Dung Tran
            % date: 4/3/2019
            % update: 11/24/2020
            
            switch nargin
                case 2
                    I = varargin{1};
                    gamma = varargin{2};
                    dis_opt = [];
                    lp_solver = 'linprog';
                case 3
                    I = varargin{1};
                    gamma = varargin{2};
                    dis_opt = varargin{3};
                    lp_solver = 'linprog';
                case 4
                    I = varargin{1};
                    gamma = varargin{2};
                    dis_opt = varargin{3};
                    lp_solver = varargin{4};
                otherwise
                    error('Invalid number of input arguments, should be 2, 3 or 4');
            end
            
            if ~isa(I, 'Star')
                error('Input is not a star');
            end

            if isempty(I)
                S = [];
            else
                [lb, ub] = I.estimateRanges;
                
                map = find(ub <= 0); % computation map
                V = I.V;
                V(map, :) = gamma*V(map, :);
                if ~isempty(I.Z)
                    c1 = I.Z.c;
                    c1(map) = gamma*c1(map);
                    V1 = I.Z.V;
                    V1(map, :) = gamma*V1(map, :);
                    new_Z = Zono(c1, V1);
                else
                    new_Z = [];
                end
                In = Star(V, I.C, I.d, I.predicate_lb, I.predicate_ub, new_Z);                    
                map = find(lb < 0 & ub > 0);
                m = length(map); 
                for i=1:m
                    if strcmp(dis_opt, 'display')
                        fprintf('\nPerforming approximate PosLin_%d operation using Star', map(i));
                    end
                    In = LeakyReLU.stepReachStarApprox(In, map(i), gamma, lp_solver);
                end
                S = In;
             
            end

        end
        
        % step reach approximation using star
        function S = multipleStepReachStarApprox_at_one(I, gamma, index, lb, ub)
            % @I: star set input
            % @gamma: coefficence of leakyReLU
            % @index: index of the neurons performing stepReach
            % @lb:lower bound of x[index]
            % @ub: upper bound of x[index]

            % author: Dung Tran
            % date: 11/24/2020
                       
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
            
            % case 0: keep the old constraints on the old predicate variables
            n = I.nVar; % number of old predicate variables
            C0 = [I.C zeros(size(I.C, 1), m)];
            d0 = I.d; 
            %case 1: y(index) >= gamma*x(index)           
            C1 = [gamma*I.V(index, 2:n+1) -V2(index, 1:m)];
            d1 = -gamma*I.V(index, 1);
            %case 2: y(index) >= x(index)
            C2 = [I.V(index, 2:n+1) -V2(index, 1:m)];
            d2 = -I.V(index, 1);
            %case 3: y(index) <= [(ub-gamma*lb)/(ub-lb)](x-lb) + gamma*lb
            a = (ub - gamma*lb)./(ub - lb); % divide element-wise
            b = a.*lb; % multiply element-wise
            C3 = [-a.*I.V(index, 2:n+1) V2(index, 1:m)];
            d3 = a.*I.V(index, 1)-b + gamma*lb;

            % create new Star
            new_C = [C0; C1; C2; C3];
            new_d = [d0; d1; d2; d3];
            new_pred_lb = [I.predicate_lb; gamma*lb];
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
            % date: 11/24/2020         
            
            switch nargin
                case 2
                    I = varargin{1};
                    gamma = varargin{2};
                    option = 'single';
                    dis_opt = [];
                    lp_solver = 'linprog';
                case 3
                    I = varargin{1};
                    gamma = varargin{2};
                    option = varargin{3};
                    dis_opt = [];
                    lp_solver = 'linprog';
                case 4
                    I = varargin{1};
                    gamma = varargin{2};
                    option = varargin{3};
                    dis_opt = varargin{4};
                    lp_solver = 'linprog';
                case 5
                    I = varargin{1};
                    gamma = varargin{2};
                    option = varargin{3};
                    dis_opt = varargin{4};
                    lp_solver = varargin{5};
                otherwise
                    error('Invalid number of input arguments, should be 2, 3, 4, or 5');
            end

            if ~isa(I, 'Star')
                error('Input is not a star');
            end

            if isempty(I)
                S = [];
            else
                [lb, ub] = I.estimateRanges;
               
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
                In = I.scaleRow(map11, gamma); % scale the element having ub <= 0
                if strcmp(dis_opt, 'display')
                    fprintf('\n(%d+%d =%d)/%d neurons have ub <= 0', length(map1), length(map3), length(map11), length(ub));
                end

                % find all indexes that have lb < 0 & ub > 0, then apply the over-approximation rule for leakyReLU
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
                S = LeakyReLU.multipleStepReachStarApprox_at_one(In, gamma, map8, lb1, ub1); % one-shot approximation
            end

        end        
        
    end
    
    methods(Static) % reachability analysis using relax-star method

        % a relaxed star-approx method using distance heuristics
        function S = reach_relaxed_star_dis(varargin)
            % @I: star input set
            % @relaxFactor: a relaxFactor
            % @S: star output set

            % author: Dung Tran
            % date: 11/24/2020
            
            switch nargin
                case 2
                    I = varargin{1};
                    gamma = varargin{2};
                    relaxFactor = 0;
                    option = 'single';
                    dis_opt = [];
                    lp_solver = 'linprog';
                case 3
                    I = varargin{1};
                    gamma = varargin{2};
                    relaxFactor = varargin{3};
                    option = 'single';
                    dis_opt = [];
                    lp_solver = 'linprog';
                case 4
                    I = varargin{1};
                    gamma = varargin{2};
                    relaxFactor = varargin{3};
                    option = varargin{4};
                    dis_opt = [];
                    lp_solver = 'linprog';
                case 5
                    I = varargin{1};
                    gamma = varargin{2};
                    relaxFactor = varargin{3};
                    option = varargin{4};
                    dis_opt = varargin{5};
                    lp_solver = 'linprog';
                case 6
                    I = varargin{1};
                    gamma = varargin{2};
                    relaxFactor = varargin{3};
                    option = varargin{4};
                    dis_opt = varargin{5};
                    lp_solver = varargin{6};
                otherwise
                    error('Invalid number of input arguments, should be 2, 3, 4, 5 or 6');
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
                    In = I.scaleRow(map11, gamma); % reset to zero at the element having ub <= 0
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
                    S = LeakyReLU.multipleStepReachStarApprox_at_one(In, gamma, map9, lb3, ub3); % one-shot approximation
                end
            end

        end
        
        % a relaxed star-approx method using area heuristic
        function S = reach_relaxed_star_area(varargin)
            % @I: star input set
            % @relaxFactor: a relaxFactor
            % @S: star output set
            % % optimize ranges of neurons that have largest estimated areas

            % author: Dung Tran
            % date: 11/24/2020
            
            switch nargin
                case 2
                    I = varargin{1};
                    gamma = varargin{2};
                    relaxFactor = 0;
                    option = 'single';
                    dis_opt = [];
                    lp_solver = 'linprog';
                case 3
                    I = varargin{1};
                    gamma = varargin{2};
                    relaxFactor = varargin{3};
                    option = 'single';
                    dis_opt = [];
                    lp_solver = 'linprog';
                case 4
                    I = varargin{1};
                    gamma = varargin{2};
                    relaxFactor = varargin{3};
                    option = varargin{4};
                    dis_opt = [];
                    lp_solver = 'linprog';
                case 5
                    I = varargin{1};
                    gamma = varargin{2};
                    relaxFactor = varargin{3};
                    option = varargin{4};
                    dis_opt = varargin{5};
                    lp_solver = 'linprog';
                case 6
                    I = varargin{1};
                    gamma = varargin{2};
                    relaxFactor = varargin{3};
                    option = varargin{4};
                    dis_opt = varargin{5};
                    lp_solver = varargin{6};
                otherwise
                    error('Invalid number of input arguments, should be 2, 3, 4, 5 or 6');
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
                
                areas = 0.5*(abs(ub(map2)).*abs(lb(map2))); % estimated areas of triangle overapproximation at all neurons
                [~,midx] = sort(areas, 'descend');
                map21 = map2(midx(1:n1)); % neurons with optimized ranged
                map22 = map2(midx(n1+1:length(map2))); % neurons without optimized ranges
                lb1 = lb(map22);
                ub1 = ub(map22); 

                if strcmp(dis_opt, 'display')
                    fprintf('\nOptimize upper bounds of %d neurons: ', length(map21));
                end
                xmax = I.getMaxs(map21, option, dis_opt, lp_solver); 
                map3 = find(xmax <= 0);

                n = length(map3);
                map4 = zeros(n,1);
                for i=1:n
                    map4(i) = map21(map3(i));
                end
                map11 = [map1; map4];
                In = I.scaleRow(map11, gamma); % reset to zero at the element having ub <= 0
                
                % find all indexes that have lb < 0 & ub > 0, then
                % apply the over-approximation rule for ReLU
                
                map5 = find(xmax > 0);
                map6 = map21(map5); % all indexes having ub > 0
                xmax1 = xmax(map5); % upper bound of all neurons having ub > 0
                
                if strcmp(dis_opt, 'display')
                    fprintf('\nOptimize lower bounds of %d neurons: ', length(map6));
                end
                xmin = I.getMins(map6, option, dis_opt, lp_solver); 
                map7 = find(xmin < 0); 
                map8 = map6(map7); % all indexes having lb < 0 & ub > 0
                lb2 = xmin(map7);  % lower bound of all indexes having lb < 0 & ub > 0
                ub2 = xmax1(map7); % upper bound of all neurons having lb < 0 & ub > 0

                map9 = [map22; map8];
                lb3 = [lb1; lb2];
                ub3 = [ub1; ub2];
                if strcmp(dis_opt, 'display')
                    fprintf('\n%d/%d neurons have lb < 0 & ub > 0', length(map9), length(ub));
                    fprintf('\nConstruct new star set, %d new predicate variables are introduced', length(map9));
                end
                S = LeakyReLU.multipleStepReachStarApprox_at_one(In, gamma, map9, lb3, ub3); % one-shot approximation
               
            end

        end
        
        % a relaxed star-approx method using lower bound and upper bound heuristic
        function S = reach_relaxed_star_lb_ub(varargin)
            % @I: star input set
            % @relaxFactor: a relaxFactor
            % @S: star output set
            % optimize ranges of neurons that have largest lower bounds and upper bounds

            % author: Dung Tran
            % date: 11/24/2020
            
            switch nargin
                case 2
                    I = varargin{1};
                    gamma = varargin{2};
                    relaxFactor = 0;
                    option = 'single';
                    dis_opt = [];
                    lp_solver = 'linprog';
                case 3
                    I = varargin{1};
                    gamma = varargin{2};
                    relaxFactor = varargin{3};
                    option = 'single';
                    dis_opt = [];
                    lp_solver = 'linprog';
                case 4
                    I = varargin{1};
                    gamma = varargin{2};
                    relaxFactor = varargin{3};
                    option = varargin{4};
                    dis_opt = [];
                    lp_solver = 'linprog';
                case 5
                    I = varargin{1};
                    gamma = varargin{2};
                    relaxFactor = varargin{3};
                    option = varargin{4};
                    dis_opt = varargin{5};
                    lp_solver = 'linprog';
                case 6
                    I = varargin{1};
                    gamma = varargin{2};
                    relaxFactor = varargin{3};
                    option = varargin{4};
                    dis_opt = varargin{5};
                    lp_solver = varargin{6};
                otherwise
                    error('Invalid number of input arguments, should be 2, 3, 4, 5 or 6');
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

                N = length(ub(map2));
                lu = [ub(map2); abs(lb(map2))];
                [~,midx] = sort(lu, 'descend');
                midx1 = midx(1:2*n1); % neurons with optimized ranges
                ub_idx = midx1(midx1 <= N); % neurons having upperbound optimized
                lb_idx = midx1(midx1 > N) - N;  % neurons having lowerbound optimized
                map21 = map2(ub_idx(:)); 
                map22 = map2(lb_idx(:)); 

                if strcmp(dis_opt, 'display')
                    fprintf('\nOptimize %d upper bounds of %d neurons: ', length(map21), length(map2));
                end

                if ~isempty(map21)
                    xmax = I.getMaxs(map21, option, dis_opt, lp_solver); 
                    map3 = find(xmax <= 0);
                    if strcmp(dis_opt, 'display')
                        fprintf('\n%d neurons (in %d neurons) with ub <= 0 are found by optimizing ranges', length(map3), length(map21));
                    end
                    n = length(map3);
                    map4 = zeros(n,1);
                    for i=1:n
                        map4(i) = map21(map3(i));
                    end

                    map5 = find(xmax > 0);
                    map6 = map21(map5(:));
                    map11 = [map1; map4];
                else
                    map11 = map1;
                    map5 = [];
                    map6 = [];
                    map4 = [];
                end
                In = I.scaleRow(map11, gamma); % reset to zero at the element having ub <= 0

                if strcmp(dis_opt, 'display') && ~isempty(map21)
                    fprintf('\n(%d+%d =%d)/%d neurons have ub <= 0', length(map1), length(map3), length(map11), length(ub));
                end
                
                if ~isempty(map4) 
                    map23 = setdiff(map22,map4, 'stable');
                else
                    map23 = map22;
                end

                % find all indexes that have lb < 0 & ub > 0, then
                % apply the over-approximation rule for ReLU
                if strcmp(dis_opt, 'display')
                    fprintf('\nOptimize %d lower bounds of %d neurons: ', length(map23), length(map2));
                end
                if ~isempty(map23)
                    xmin = I.getMins(map23, option, dis_opt, lp_solver); 
                    map7 = find(xmin < 0);
                    map8 = map23(map7(:));
                    map9 = find(xmin >= 0);
                    map10 = map23(map9(:));
                else
                    map8 = [];
                    map10 = [];
                end

                if ~isempty(map4)
                    map24 = setdiff(map2,map4, 'stable');
                else
                    map24 = map2;
                end

                if ~isempty(map10)
                    map24 = setdiff(map24,map10, 'stable');
                end

                if ~isempty(map6)
                    ub(map6) = xmax(map5);
                end
                if ~isempty(map8)
                    lb(map8) = xmin(map7);
                end

                ub1 = ub(map24);
                lb1 = lb(map24);

                if strcmp(dis_opt, 'display')
                    fprintf('\nConstruct new star set, %d new predicate variables are introduced', length(map24));
                end
                S = LeakyReLU.multipleStepReachStarApprox_at_one(In, gamma, map24, lb1, ub1); % one-shot approximation

            end

        end

        % a relaxed star-approx method using upper bound heuristic
        function S = reach_relaxed_star_ub(varargin)
            % @I: star input set
            % @relaxFactor: a relaxFactor
            % @S: star output set
            % optimize ranges of neurons that have largest lower bounds and upper bounds

            % author: Dung Tran
            % date: 11/24/2020
            
            switch nargin
                case 2
                    I = varargin{1};
                    gamma = varargin{2};
                    relaxFactor = 0;
                    option = 'single';
                    dis_opt = [];
                    lp_solver = 'linprog';
                case 3
                    I = varargin{1};
                    gamma = varargin{2};
                    relaxFactor = varargin{3};
                    option = 'single';
                    dis_opt = [];
                    lp_solver = 'linprog';
                case 4
                    I = varargin{1};
                    gamma = varargin{2};
                    relaxFactor = varargin{3};
                    option = varargin{4};
                    dis_opt = [];
                    lp_solver = 'linprog';
                case 5
                    I = varargin{1};
                    gamma = varargin{2};
                    relaxFactor = varargin{3};
                    option = varargin{4};
                    dis_opt = varargin{5};
                    lp_solver = 'linprog';
                case 6
                    I = varargin{1};
                    gamma = varargin{2};
                    relaxFactor = varargin{3};
                    option = varargin{4};
                    dis_opt = varargin{5};
                    lp_solver = varargin{6};
                otherwise
                    error('Invalid number of input arguments, should be 2, 3, 4, 5 or 6');
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
        
                N = length(ub(map2));
                [~,midx_u] = sort(ub(map2), 'descend');
                [~, midx_l] = sort(abs(lb(map2)), 'descend');
                if 2*n1 <= N
                    ub_idx = midx_u(1:2*n1); 
                    map21 = map2(ub_idx);
                    map22 = [];
                else
                    ub_idx = midx_u(1:N);
                    map21 = map2(ub_idx);
                    lb_idx = midx_l(1:2*n1 - N);
                    map22 = map2(lb_idx);                 
                end              
                
                if strcmp(dis_opt, 'display')
                    fprintf('\nOptimize %d upper bounds of %d neurons: ', length(map21), length(map2));
                end

                if ~isempty(map21)
                    xmax = I.getMaxs(map21, option, dis_opt, lp_solver); 
                    map3 = find(xmax <= 0);
                    if strcmp(dis_opt, 'display')
                        fprintf('\n%d neurons (in %d neurons) with ub <= 0 are found by optimizing ranges', length(map3), length(map21));
                    end
                    n = length(map3);
                    map4 = zeros(n,1);
                    for i=1:n
                        map4(i) = map21(map3(i));
                    end
                    map5 = find(xmax > 0);
                    map6 = map21(map5(:));
                    map11 = [map1; map4];                    
                else
                    map11 = map1;
                    map5 = [];
                    map6 = [];
                    map4 = [];
                end
                In = I.scaleRow(map11, gamma); % reset to zero at the element having ub <= 0

                if strcmp(dis_opt, 'display') && ~isempty(map21)
                    fprintf('\n(%d+%d =%d)/%d neurons have ub <= 0', length(map1), length(map3), length(map11), length(ub));
                end

                if ~isempty(map22)
                    if ~isempty(map4)
                        map23 = setdiff(map22,map4, 'stable');
                    else
                        map23 = map22;
                    end
                else
                    map23 = map22;
                end
                

                % find all indexes that have lb < 0 & ub > 0, then
                % apply the over-approximation rule for ReLU
                if strcmp(dis_opt, 'display')
                    fprintf('\nOptimize %d lower bounds of %d neurons: ', length(map23), length(map2));
                end
                if ~isempty(map23)
                    xmin = I.getMins(map23, option, dis_opt, lp_solver); 
                    map7 = find(xmin < 0);
                    map8 = map23(map7(:));
                    map9 = find(xmin >= 0);
                    map10 = map23(map9(:));
                else
                    map8 = [];
                    map10 = [];
                end

                if ~isempty(map4)
                    map24 = setdiff(map2, map4, 'stable');
                else
                    map24 = map2;
                end

                if ~isempty(map10)
                    map24 = setdiff(map24, map10, 'stable');
                end

                if ~isempty(map6)
                    ub(map6) = xmax(map5);
                end
                if ~isempty(map8)
                    lb(map8) = xmin(map7);
                end

                ub1 = ub(map24);
                lb1 = lb(map24);

                if strcmp(dis_opt, 'display')
                    fprintf('\nConstruct new star set, %d new predicate variables are introduced', length(map24));
                end
                S = LeakyReLU.multipleStepReachStarApprox_at_one(In, gamma, map24, lb1, ub1); % one-shot approximation

            end

        end
        
        % a relaxed star-approx method using random heuristic
        function S = reach_relaxed_star_random(varargin)
            % @I: star input set
            % @relaxFactor: a relaxFactor
            % @S: star output set
            % optimize ranges of randomly selected neurons

            % author: Dung Tran
            % date: 11/24/2020
            
            switch nargin
                case 2
                    I = varargin{1};
                    gamma = varargin{2};
                    relaxFactor = 0;
                    option = 'single';
                    dis_opt = [];
                    lp_solver = 'linprog';
                case 3
                    I = varargin{1};
                    gamma = varargin{2};
                    relaxFactor = varargin{3};
                    option = 'single';
                    dis_opt = [];
                    lp_solver = 'linprog';
                case 4
                    I = varargin{1};
                    gamma = varargin{2};
                    relaxFactor = varargin{3};
                    option = varargin{4};
                    dis_opt = [];
                    lp_solver = 'linprog';
                case 5
                    I = varargin{1};
                    gamma = varargin{2};
                    relaxFactor = varargin{3};
                    option = varargin{4};
                    dis_opt = varargin{5};
                    lp_solver = 'linprog';
                case 6
                    I = varargin{1};
                    gamma = varargin{2};
                    relaxFactor = varargin{3};
                    option = varargin{4};
                    dis_opt = varargin{5};
                    lp_solver = varargin{6};
                otherwise
                    error('Invalid number of input arguments, should be 2, 3, 4, 5 or 6');
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

                midx = randperm(length(map2), n1);
                midx = midx';

                map21 = map2(midx(1:n1)); % neurons with optimized ranged
                map22 = setdiff(map2, map21, 'stable');
                lb1 = lb(map22);
                ub1 = ub(map22); 

                if strcmp(dis_opt, 'display')
                    fprintf('\nOptimize upper bounds of %d neurons: ', length(map21));
                end
                xmax = I.getMaxs(map21, option, dis_opt, lp_solver); 
                map3 = find(xmax <= 0);

                n = length(map3);
                map4 = zeros(n,1);
                for i=1:n
                    map4(i) = map21(map3(i));
                end
                map11 = [map1; map4];
                In = I.scaleRow(map11, gamma); % reset to zero at the element having ub <= 0
                
                % find all indexes that have lb < 0 & ub > 0, then
                % apply the over-approximation rule for ReLU
                
                map5 = find(xmax > 0);
                map6 = map21(map5); % all indexes having ub > 0
                xmax1 = xmax(map5); % upper bound of all neurons having ub > 0
                
                if strcmp(dis_opt, 'display')
                    fprintf('\nOptimize lower bounds of %d neurons: ', length(map6));
                end
                xmin = I.getMins(map6, option, dis_opt, lp_solver); 
                map7 = find(xmin < 0); 
                map8 = map6(map7); % all indexes having lb < 0 & ub > 0
                lb2 = xmin(map7);  % lower bound of all indexes having lb < 0 & ub > 0
                ub2 = xmax1(map7); % upper bound of all neurons having lb < 0 & ub > 0
                
                map9 = [map22; map8];
                lb3 = [lb1; lb2];
                ub3 = [ub1; ub2];
                if strcmp(dis_opt, 'display')
                    fprintf('\n%d/%d neurons have lb < 0 & ub > 0', length(map9), length(ub));
                    fprintf('\nConstruct new star set, %d new predicate variables are introduced', length(map9));
                end
                S = LeakyReLU.multipleStepReachStarApprox_at_one(In, gamma, map9, lb3, ub3); % one-shot approximation
               
            end

        end
        
        % a relaxed star-approx method using static heuristic
        function S = reach_relaxed_star_static(varargin)
            % @I: star input set
            % @relaxFactor: a relaxFactor
            % @S: star output set
            % optimize ranges of the first n neurons

            % author: Dung Tran
            % date: 11/24/2020
            
            switch nargin
                case 2
                    I = varargin{1};
                    gamma = varargin{2};
                    relaxFactor = 0;
                    option = 'single';
                    dis_opt = [];
                    lp_solver = 'linprog';
                case 3
                    I = varargin{1};
                    gamma = varargin{2};
                    relaxFactor = varargin{3};
                    option = 'single';
                    dis_opt = [];
                    lp_solver = 'linprog';
                case 4
                    I = varargin{1};
                    gamma = varargin{2};
                    relaxFactor = varargin{3};
                    option = varargin{4};
                    dis_opt = [];
                    lp_solver = 'linprog';
                case 5
                    I = varargin{1};
                    gamma = varargin{2};
                    relaxFactor = varargin{3};
                    option = varargin{4};
                    dis_opt = varargin{5};
                    lp_solver = 'linprog';
                case 6
                    I = varargin{1};
                    gamma = varargin{2};
                    relaxFactor = varargin{3};
                    option = varargin{4};
                    dis_opt = varargin{5};
                    lp_solver = varargin{6};
                otherwise
                    error('Invalid number of input arguments, should be 2, 3, 4, 5 or 6');
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
                                
                map21 = map2(1:n1); % neurons with optimized ranged
                map22 = map2(n1+1:length(map2));
                lb1 = lb(map22);
                ub1 = ub(map22); 

                if strcmp(dis_opt, 'display')
                    fprintf('\nOptimize upper bounds of %d neurons: ', length(map21));
                end
                xmax = I.getMaxs(map21, option, dis_opt, lp_solver); 
                map3 = find(xmax <= 0);

                n = length(map3);
                map4 = zeros(n,1);
                for i=1:n
                    map4(i) = map21(map3(i));
                end
                map11 = [map1; map4];
                In = I.scaleRow(map11, gamma); % reset to zero at the element having ub <= 0
                
                % find all indexes that have lb < 0 & ub > 0, then
                % apply the over-approximation rule for ReLU
                
                map5 = find(xmax > 0);
                map6 = map21(map5); % all indexes having ub > 0
                xmax1 = xmax(map5); % upper bound of all neurons having ub > 0
                
                if strcmp(dis_opt, 'display')
                    fprintf('\nOptimize lower bounds of %d neurons: ', length(map6));
                end
                xmin = I.getMins(map6, option, dis_opt, lp_solver); 
                map7 = find(xmin < 0); 
                map8 = map6(map7); % all indexes having lb < 0 & ub > 0
                lb2 = xmin(map7);  % lower bound of all indexes having lb < 0 & ub > 0
                ub2 = xmax1(map7); % upper bound of all neurons having lb < 0 & ub > 0

                map9 = [map22; map8];
                lb3 = [lb1; lb2];
                ub3 = [ub1; ub2];
                if strcmp(dis_opt, 'display')
                    fprintf('\n%d/%d neurons have lb < 0 & ub > 0', length(map9), length(ub));
                    fprintf('\nConstruct new star set, %d new predicate variables are introduced', length(map9));
                end
                S = LeakyReLU.multipleStepReachStarApprox_at_one(In, gamma, map9, lb3, ub3); % one-shot approximation
               
            end

        end
        
    end
    
    
    methods(Static) % reachability analysis using Polyhedron
        
        % step reach set y = leakyReLU(x)
        function R = stepReach_Polyhedron(I, gamma, index, xmin, xmax)
            % @I : input set, a polyhedron
            % @gamma: coefficience of leakyReLU
            % @i : index of current x[index] of current step
            % @xmin: min of x[index]
            % @xmax: max of x[index]
            
            % author: Dung Tran
            % date: 11/24/2020
           
            I.normalize;
            dim = I.Dim;
            if xmin >= 0
                R = I; 
            elseif xmax < 0 
                Im = eye(dim);
                Im(index, index) = gamma;
                R = Im*I;
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
                Im(index, index) = gamma;
                R1 = Im*R1;
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
                
                case 6
                    I = varargin{1};
                    gamma = varargin{2};
                    index = varargin{3};
                    xmin = varargin{4};
                    xmax = varargin{5};
                    option = varargin{6};
                
                otherwise
                    error('Invalid number of input arguments (should be 6)');
            end
               
            p = length(I);
            R = [];
            
            if isempty(option)
                
                for i=1:p
                    R =[R, LeakyReLU.stepReach_Polyhedron(I(i), gamma, index, xmin, xmax)];
                end
                
            elseif strcmp(option, 'parallel')
                
                parfor i=1:p
                    R =[R, LeakyReLU.stepReach_Polyhedron(I(i), gamma, index, xmin, xmax)];
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
                    gamma = varargin{2};
                    option = [];
                    dis_opt = [];
                case 3
                    I = varargin{1};
                    gamma = varargin{2};
                    option = varargin{3};
                    dis_opt = [];
                case 4
                    I = varargin{1};
                    gamma = varargin{2};
                    option = varargin{3};
                    dis_opt = varargin{4};
                otherwise
                    error('Invalid number of input arguments, should be 3 or 4');
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
                        In = LeakyReLU.stepReachMultipleInputs_Polyhedron(In, gamma, map(i), lb(map(i)), ub(map(i)), option);
                    end               
                    R = In;
                end
                
            else
                R = [];
            end
                 
        end
             
    end


    methods(Static) % over-approximate reachability analysis use zonotope            
            
        % over-approximate reachability analysis use zonotope
        function Z = reach_zono_approx(varargin)
            % @I: zonotope input
            % @Z: zonotope output
            
            % author: Dung Tran
            % date: 11/24/2020
            
            % reference: Fast and Effective Robustness Ceritification,
            % Gagandeep Singh, NIPS 2018
            
            switch nargin
                case 2
                    I = varargin{1};
                    gamma = varargin{2};
                    dis_opt = [];
                case 3
                    I = varargin{1};
                    gamma = varargin{2};
                    dis_opt = varargin{3};
                otherwise
                    error('Invalid number of input arguments, should be 2 or 3');
            end
            
            if ~isa(I, 'Zono')
                error('Input is not a Zonotope');
            end
                      
            [lb, ub] = I.getBounds;
            map1 = find(ub <= 0);
            c1 = I.c;
            V1 = I.V;
            c1(map1) = gamma*c1(map1);
            V1(map1, :) = gamma*V1(map1);
            map2 = find(ub > 0 & lb < 0);
            a = (ub(map2) - gamma*lb(map2))./(ub(map2) - lb(map2));
            mu = gamma*lb(map2) - a.*lb(map2);
            
            c1(map2) = a.*c1(map2) + mu;
            V1(map2, :) = a.*V1(map2, :);
            m = length(map2);
            V2 = zeros(I.dim, m);
            for i=1:m
                V2(map2(i), i) = 1;
            end
            V2(map2, :) = mu.*V2(map2, :);
            if strcmp(dis_opt, 'display')
                fprintf('\nPerforming %d approximate leakyReLU operations using Zonotope', m);
            end
            Z = Zono(c1, [V1 V2]);
                       
        end
        
    end
    
    
    methods(Static) % reachability analysis using abstract-domain
        
        % step over-approximate reachability analysis using abstract-domain
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
                
                case 5
                    I = varargin{1};
                    gamma = varargin{2};
                    index = varargin{3};
                    lb = varargin{4};
                    ub = varargin{5};
                otherwise
                    error('Invalid number of input arguments (should be 5)');
            end
            
            if ~isa(I, 'Star')
                error('Input is not a Star');
            end
                          
            if lb >= 0
                S = Star(I.V, I.C, I.d, I.predicate_lb, I.predicate_ub);
            elseif ub <= 0
                V = I.V;
                V(index, :) = gamma*V(index, :);
                S = Star(V, I.C, I.d, I.predicate_lb, I.predicate_ub);
            elseif lb < 0 && ub > 0
                
                dis1 = abs(ub*(gamma-1))/(sqrt(gamma^2 + 1));
                len1 = abs(lb)*sqrt(gamma^2 +1);
                dis2 = abs(lb*(1-gamma))/sqrt(2);
                len2 = ub*sqrt(2);
                S1 = dis1*len1; % area of the first candidate abstract-domain
                S2 = dis2*len2; % area of the second candidate abstract-domain  
                                
                n = I.nVar + 1;
                                
                % constraint 1: y[index] = leakyReLU(x[index]) >= gamma*x[index]
                C1 = [gamma*I.V(index, 2:n) -1];
                d1 = -gamma*I.V(index, 1);
                
                % constraint 2: y[index] = ReLU(x[index]) >= x[index]
                C2 = [I.V(index, 2:n) -1];
                d2 = -I.V(index, 1);
                    
                % constraint 3: y[index] <=
                % [(ub-gamma*lb)/(ub-lb)](x-lb) + gamma*lb
                a = (ub-gamma*lb)/(ub-lb);
                C3 = [-a*I.V(index, 2:n) 1];
                d3 = a*I.V(index, 1) - a*lb + gamma*lb;
                               
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
                    new_pred_lb = [I.predicate_lb; gamma*lb];
                    new_pred_ub = [I.predicate_ub; ub];
                    
                else
                    % choose the second candidate as the abstract-domain                                      
                    new_C = [C0; C2; C3];
                    new_d = [d0; d2; d3];
                    new_pred_lb = [I.predicate_lb; lb];
                    new_pred_ub = [I.predicate_ub; ub];
                                        
                end
                
                S = Star(new_V, new_C, new_d, new_pred_lb, new_pred_ub);
                
            end
                       
        end
        
        % over-approximate reachability analysis using abstract-domain
        function S = reach_abstract_domain(varargin)
            % @I: star input set
            % @S: star output set

            % author: Dung Tran
            % date: 11/24/2020
            
            switch nargin
                case 2
                    I = varargin{1};
                    gamma = varargin{2};
                    dis_opt = [];
                    lp_solver = 'linprog';
                case 3
                    I = varargin{1};
                    gamma = varargin{2};
                    dis_opt = varargin{3};
                    lp_solver = 'linprog';
                case 4
                    I = varargin{1};
                    gamma = varargin{2};
                    dis_opt = varargin{3};
                    lp_solver = varargin{4};
                otherwise
                    error('Invalid number of input arguments, should be 2, 3, or 4');
            end
            
            if ~isa(I, 'Star')
                error('Input is not a star');
            end

            if isempty(I)
                S = [];
            else    
                %[lb, ub] = I.estimateRanges;
                
                % TODO: minimize the number of LP needs to be solved
                map = 1:1:I.dim;
                lb = I.getMins(map, 'single', lp_solver);
                ub = I.getMaxs(map, 'single', lp_solver);
                if isempty(lb) || isempty(ub)
                    S = [];
                else
                    map = find(ub <= 0); % computation map
                    V = I.V;
                    V(map, :) = gamma*V(map, :); 
                    In = Star(V, I.C, I.d, I.predicate_lb, I.predicate_ub);                    
                    map = find(lb < 0 & ub > 0);
                    m = length(map); 
                    for i=1:m
                        if strcmp(dis_opt, 'display')
                            fprintf('\nPerforming approximate PosLin_%d operation using Abstract Domain', map(i));
                        end
                        In = LeakyReLU.stepReachAbstractDomain(In, gamma, map(i), lb(map(i)), ub(map(i)));
                    end
                    S = In;
             
                end
            
            end     
        
        end
        
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
                
                case 7
                    I = varargin{1};
                    gamma = varargin{2};
                    method = varargin{3};
                    option = varargin{4};
                    relaxFactor = varargin{5}; % used for aprox-star only
                    dis_opt = varargin{6}; % display option
                    lp_solver = varargin{7}; 
                case 6
                    I = varargin{1};
                    gamma = varargin{2};
                    method = varargin{3};
                    option = varargin{4};
                    relaxFactor = varargin{5}; % used for aprox-star only
                    dis_opt = varargin{6}; % display option
                    lp_solver = 'linprog';
                case 5
                    I = varargin{1};
                    gamma = varargin{2};
                    method = varargin{3};
                    option = varargin{4};
                    relaxFactor = varargin{5}; % used for aprox-star only
                    dis_opt = [];
                    lp_solver = 'linprog';
                case 4
                    I = varargin{1};
                    gamma = varargin{2};
                    method = varargin{3};
                    option = varargin{4};
                    relaxFactor = 0; % used for aprox-star only
                    dis_opt = [];
                    lp_solver = 'linprog';
                case 3
                    I = varargin{1};
                    gamma = varargin{2};
                    method = varargin{3};
                    option = [];
                    relaxFactor = 0; % used for aprox-star only
                    dis_opt = [];
                    lp_solver = 'linprog';
                case 2
                    I = varargin{1};
                    gamma = varargin{2};
                    method = 'exact-star';
                    option = [];
                    relaxFactor = 0; % used for aprox-star only
                    dis_opt = [];
                    lp_solver = 'linprog';
                otherwise
                    error('Invalid number of input arguments (should be 2,3,4,5,6 or 7)');
            end
            
            if strcmp(method, 'exact-star') % exact analysis using star
                R = LeakyReLU.reach_star_exact(I, gamma, option, dis_opt, lp_solver);
            elseif strcmp(method, 'exact-polyhedron') % exact analysis using polyhedron
                R = LeakyReLU.reach_polyhedron_exact(I, gamma, option, dis_opt);
            elseif strcmp(method, 'approx-star')  % over-approximate analysis using star
                R = LeakyReLU.reach_star_approx2(I, gamma, option, dis_opt, lp_solver);
            elseif strcmp(method, 'relax-star-dis')
                R = LeakyReLU.reach_relaxed_star_dis(I, gamma, relaxFactor, option, dis_opt, lp_solver);
            elseif strcmp(method, 'relax-star-lb-ub')
                R = LeakyReLU.reach_relaxed_star_lb_ub(I, gamma, relaxFactor, option, dis_opt, lp_solver);
            elseif strcmp(method, 'relax-star-area')
                R = LeakyReLU.reach_relaxed_star_area(I, gamma, relaxFactor, option, dis_opt, lp_solver);
            elseif strcmp(method, 'relax-star-ub')
                R = LeakyReLU.reach_relaxed_star_ub(I, gamma, relaxFactor, option, dis_opt, lp_solver);
            elseif strcmp(method, 'relax-star-random')
                R = LeakyReLU.reach_relaxed_star_random(I, gamma, relaxFactor, option, dis_opt, lp_solver);
            elseif strcmp(method, 'relax-star-static')
                R = LeakyReLU.reach_relaxed_star_static(I, gamma, relaxFactor, option, dis_opt, lp_solver);
            elseif strcmp(method, 'approx-zono')  % over-approximate analysis using zonotope 
                R = LeakyReLU.reach_zono_approx(I, gamma, dis_opt);
            elseif strcmp(method, 'abs-dom')  % over-approximate analysis using abstract-domain
                R = LeakyReLU.reach_abstract_domain(I, gamma, dis_opt, lp_solver);
            elseif strcmp(method, 'exact-face-latice') % exact analysis using face-latice
                fprintf('\nNNV have not yet support Exact Face-Latice Method');
            elseif strcmp(method, 'approx-face-latice') % over-approximate analysis using face-latice
                fprintf('\nNNV have not yet support Approximate Face-Latice Method');
            end
                            
        end
        
    end
    
end
