classdef Sign
    % Sign Class contains method for reachability analysis for Layer with sign activation function 
    % Its name is consistent with the name of ReLU acivation function in
    % matlab called 'sign'
    
    % author: Mykhailo Ivashchenko
    
    properties(Constant)                   
        sign_bounds = containers.Map({'polar_zero_to_pos_one', 'nonnegative_zero_to_pos_one'}, {[-1, 1], [0, 1]});
    end
    
    methods(Static) % evaluate method and reachability analysis with stars
       
        function sign_lb = get_sign_lb(map)
            sign_lb = Sign.sign_bounds(map);
            sign_lb = sign_lb(1);
        end
        
        function sign_ub = get_sign_ub(map)
            sign_ub = Sign.sign_bounds(map);
            sign_ub = sign_ub(2);
        end
        
        % evaluation
        function y = evaluate(x, mode)
            y = sign(x);
            
            if strcmp(mode, "polar_zero_to_pos_one")
                y = y + (y == 0);
            elseif strcmp(mode, "nonnegative_zero_to_pos_one")
                y = y + (y == 0);
                y = y + (y == -1);
            end
        end
        
        % stepReach method, compute reachable set for a single step
        function S = stepReach(varargin)
            % @I: single star set input
            % @index: index of the neuron performing stepSign
            % @xmin: minimum of x[index]
            % @xmax: maximum of x[index]
            % @S: star output set
            
            % author: Mykhailo Ivashchenko
            
            switch nargin
                case 2 
                    I = varargin{1};
                    index = varargin{2};
                    lp_solver = 'linprog';
                case 3
                    I = varargin{1};
                    index = varargin{2};
                    lp_solver = varargin{3};
                case 4
                    I = varargin{1};
                    index = varargin{2};
                    lp_solver = varargin{3};
                    mode = varargin{4};
                otherwise
                    error('Invalid number of input arguments, should be 2 or 3');
            end
            
            if ~isa(I, 'Star')
                error('Input is not a star set');
            end
            
            xmin = I.getMin(index, lp_solver);
            xmax = I.getMax(index, lp_solver);
            
            sign_lb = Sign.get_sign_lb(mode);
            sign_ub = Sign.get_sign_ub(mode);
            sign = sign_lb;
            %sign = -1;
            
            %if xmin == xmax && ((xmin == -1) || (xmin == 1))
            if xmin == xmax && ((xmin == sign_lb) || (xmin == sign_ub))
                S = I;
                return;
            %elseif xmin >= 0
            elseif xmin >= sign_lb
                sign = sign + sign_ub - sign_lb;
            end
            
            % 1. Sign(I) = -1 if sign == -1; next if statement is executed to
            % obtain Sign(I) = 1
            % 2. Sign(I) = 1 if sign == 1; next if statement is skipped
            new_V = I.V;
            new_V(index, 1) = sign;

            new_predicate_lb = I.predicate_lb;
            new_predicate_ub = I.predicate_ub;
                        
            new_predicate_lb(index) = sign;
            new_predicate_ub(index) = sign;
            new_V(index, 2:end) = zeros(1, size(new_V, 2) - 1); 
            
            new_C = I.C;
            new_d = I.d;            
            
            
            % update outer-zono
            if ~isempty(I.Z)
                c1 = I.Z.c;
                c1(index) = sign;
                V1 = I.Z.V;
                
                V1(index, 2:end) = zeros(1, size(V1, 2) - 1);
                
                new_Z = Zono(c1, V1);
            else
                new_Z = [];
            end
            S1 = Star(new_V, new_C, new_d, new_predicate_lb, new_predicate_ub, new_Z); 
            
            S = S1;
            
            if sign ~= sign_ub && xmax > 0 % S1 = I && x[id] >= 0
                sign = 1;
                new_V(index, 1) = sign;

                new_predicate_lb(index) = sign;
                new_predicate_ub(index) = sign;

                
                % update outer-zono
                if ~isempty(I.Z)
                    new_Z.c(index) = sign;
                end

                S1 = Star(new_V, new_C, new_d, new_predicate_lb, new_predicate_ub, new_Z); 

                S = [S S1];
            end

        end
        
        
        % stepReach with multiple inputs
        function S = stepReachMultipleInputs(varargin)
            % @I: an array of stars
            % @index: index where stepReach is performed
            % @option: = 'parallel' use parallel computing
            %          = not declare -> don't use parallel computing
            
            % author: Mykhailo Ivashchenko
            
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
                case 5
                    I = varargin{1};
                    index = varargin{2};
                    option = varargin{3};
                    lp_solver = varargin{4};
                    mode = varargin{5};
                otherwise
                    error('Invalid number of input arguments');
            end
            
            p = length(I);
            S = [];
            
            if isempty(option)
                
                for i=1:p
                    S =[S, Sign.stepReach(I(i), index, lp_solver, mode)];
                end
                
            elseif strcmp(option, 'parallel')
                
                parfor i=1:p
                    S =[S, Sign.stepReach(I(i), index, lp_solver, mode)];
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
            
            % author: Mykhailo Ivashchenko
            
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
                case 5
                    I = varargin{1};
                    option = varargin{2};
                    dis_opt = varargin{3};
                    lp_solver = varargin{4};
                    mode = varargin{5};
                otherwise
                    error('Invalid number of input arguments, should be 2, 3 or 4');
            end
            
             if ~isempty(I)
                if isa(I, 'ImageStar')
                    Converted_star = I.toStar();
                    In = Star(Converted_star.V, Converted_star.C, Converted_star.d, Converted_star.predicate_lb, Converted_star.predicate_ub, Converted_star.Z);
                else
                    In = Star(I.V, I.C, I.d, I.predicate_lb, I.predicate_ub, I.Z);
                end

                if(isempty(In.predicate_lb) || isempty(In.predicate_ub))
                    [In.predicate_lb, In.predicate_ub] = In.getPredicateBounds();
                end
                
                [lb, ub] = In.estimateRanges;

                if isempty(lb) || isempty(ub)
                    S = [];
                else
                    m = length(lb);                    
                    for i=1:m
                        if strcmp(dis_opt, 'display')

                            fprintf('\nPerforming exact Sign_%d operation using Star', i);
                        end
                        if i == 67
                            disp(i);
                        end
                        In = Sign.stepReachMultipleInputs(In, i, option, lp_solver, mode);
                    end               
                    S = In;
                end
                
            else
                S = [];
            end
            
        end
                
        
        % exact reachability analysis using star
        function S = reach_star_exact_multipleInputs(varargin)
            % @In: star input sets
            % @option: = 'parallel' using parallel computing
            %          = '[]'    do not use parallel computing
            
            % author: Mykhailo Ivashchenko
            
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
                case 5
                    In = varargin{1};
                    option = varargin{2};
                    dis_opt = varargin{3};
                    lp_solver = varargin{4};
                    mode = varargin{5};
                otherwise
                    error('Invalid number of input arguments, should be 2, 3 or 4');
            end
            
            
             n = length(In);
             S = [];
             if strcmp(option, 'parallel')
                 parfor i=1:n
                     S = [S Sign.reach_star_exact(In(i), [], dis_opt, lp_solver)];
                 end
             elseif isempty(option) || strcmp(option, 'single')
                 for i=1:n
                     S = [S Sign.reach_star_exact(In(i), [], dis_opt, lp_solver, mode)];
                 end
             else
                 error('Unknown computation option');
             end
        end
              
        
        
        % step reach approximation using star
        function [S, shutdownIds] = stepReachStarApprox(I, index, shutdownIds)
            % @I: star set input
            % @index: index of the neuron performing stepReach
            % @S: star output

            % author: Mykhailo Ivashchenko
                       
            
            if ~isa(I, 'Star')
                error('Input is not a star');
            end
            
            [lb, ub] = I.estimateRange(index);

            new_predicate_lb = I.predicate_lb;
            new_predicate_ub = I.predicate_ub;
            
            n = I.nVar;
            
            new_V = I.V;
            
            C0 = I.C;
            d0 = I.d;
              
            if index > I.nVar
                n = I.nVar + 1;
                
                new_V = [new_V zeros(size(new_V, 1), 1)];
                C0 = [C0 0];
                
                new_predicate_lb = [new_predicate_lb; 0];
                new_predicate_ub = [new_predicate_ub; 0];
            end
            
            if lb >= 0
                shutdownIds = [shutdownIds index];

                new_V(index, 2:end) = zeros(1, size(new_V, 2) - 1);
                new_V(index, 1) = 1;

                new_predicate_lb(index) = 1;
                new_predicate_ub(index) = 1;
                
                if ~isempty(I.Z)
                    c1 = I.Z.c;
                    c1(index) = Sign.evaluate(c1(index));

                    V1 = new_V(:, 2:end);

                    V1(index, :) = zeros(1, size(V1, 2));

                    new_Z = Zono(c1, V1);
                else
                    new_Z = [];
                end
            elseif ub < 0   
                shutdownIds = [shutdownIds index];

                new_V(index, 2:end) = zeros(1, size(new_V, 2) - 1);
                new_V(index, 1) = -1;

                new_predicate_lb(index) = -1;
                new_predicate_ub(index) = -1;
                
                if ~isempty(I.Z)
                    c1 = I.Z.c;
                    c1(index) = Sign.evaluate(c1(index));

                    V1 = new_V(:, 2:end);

                    V1(index, :) = zeros(1, size(V1, 2));

                    new_Z = Zono(c1, V1);
                else
                    new_Z = [];
                end
            else                
                new_V(index, 2 : end) = zeros(1, size(new_V, 2) - 1);
                
                new_V(index, index + 1) = 1;
                
                new_V(index, 1) = Sign.evaluate(new_V(index, 1));

                new_predicate_lb(index) = -1;
                new_predicate_ub(index) = 1;
                
                if ~isempty(I.Z)
                    c1 = I.Z.c;
                    c1(index) = Sign.evaluate(c1(index));

                    V1 = new_V(:, 2:end);

                    V1(index, :) = zeros(1, size(V1, 2));
                    V1(index, index) = 1;

                    new_Z = Zono(c1, V1);
                else
                    new_Z = [];
                end
            end
            S = Star(new_V, I.C, I.d, new_predicate_lb, new_predicate_ub, new_Z);
        end
      
        function S = stepReachApproxBase(In)
            % @In: star set input
            % @index: index of the neuron performing stepReach
            % @S: star output

            % author: Mykhailo Ivashchenko
            
            [lb, ub] = In.estimateRanges();

            if isempty(lb) || isempty(ub)
                S = [];
            else
                current_nVar = In.nVar;

                new_V = In.V;
                new_C = In.C;
                new_d = In.d;
                new_predicate_lb = In.predicate_lb;
                new_predicate_ub = In.predicate_ub;
                new_Z = In.Z;

                for i = 1:size(lb, 1)

                    if (ub(i) < 0)
                        new_V(i, 2:end) = zeros(1, size(new_V, 2) - 1);

                        new_V(i,1) = -1;
                    elseif(lb(i) >= 0) 
                        new_V(i, 2:end) = zeros(1, size(new_V, 2) - 1);

                        new_V(i,1) = 1;
                    else
                        current_nVar = current_nVar + 1;

                        new_C = [new_C zeros(size(new_C, 1), 1)];

                        new_constr = zeros(2, current_nVar);
                        new_vec = zeros(2, 1);
                        
                        new_constr(1, current_nVar) = -1; % x >= -1
                        new_constr(2, current_nVar) = 1; % x < 1

                        new_vec(1) = 1;
                        new_vec(2) = 1;

                        new_predicate_lb = [new_predicate_lb; -1];
                        new_predicate_ub = [new_predicate_ub; 1];

                        new_V = [new_V zeros(size(In.V, 1), 1)];

                        new_V(i,1) = Sign.evaluate(new_V(i,1));

                        new_C = [new_C; new_constr];
                        new_d = [new_d; new_vec];
                        
                        new_V(i, 2:end) = zeros(1, size(new_V, 2) - 1);
                        new_V(i, current_nVar + 1) = 1;
                    end
                end

                if ~isempty(In.Z)
                    new_Z.c = new_V(:, 1);

                    new_Z.V = new_V(:, 2:end);
                end

                S = Star(new_V, new_C, new_d, new_predicate_lb, new_predicate_ub, new_Z);
            end
        end
      
        function S = stepReachApprox_v1(In, mode)
            % @In: star set input
            % @index: index of the neuron performing stepReach
            % @S: star output

            % author: Mykhailo Ivashchenko
            
            [lb, ub] = In.estimateRanges();
            
            sign_lb = Sign.get_sign_lb(mode);
            sign_ub = Sign.get_sign_ub(mode);

            if isempty(lb) || isempty(ub)
                S = [];
            else
                current_nVar = In.nVar;

                l = length(lb);         

                new_V = In.V;
                new_Z = In.Z;

                neg_bool = ub < 0;
                pos_bool = lb >= 0;
                if(isa(neg_bool, "dlarray"))
                    neg_bool = extractdata(neg_bool);
                    pos_bool = extractdata(pos_bool);
                end

                neg_ids = find(neg_bool);
                pos_ids = find(pos_bool);

                ids = (1 : length(ub))';
                others = setdiff(ids, neg_ids);
                others = setdiff(others, pos_ids);
                additional = length(others);
                

                new_V = [new_V zeros(l, additional)];
                new_V(:, 2 : (current_nVar + 1)) = zeros(l, current_nVar);
                new_V(others, current_nVar + 2 : end) = eye(additional, additional);
                
                %%% new_V([neg_ids; pos_ids], 2:end) = zeros(length([neg_ids; pos_ids]), size(new_V, 2) - 1);

%                 new_V(neg_ids, 1) = -1;
%                 new_V(pos_ids, 1) = 1;
                new_V(neg_ids, 1) = sign_lb;
                new_V(pos_ids, 1) = sign_ub;
                new_V(others, 1) = Sign.evaluate(new_V(others, 1), mode);

                new_C = zeros((2 * additional), current_nVar + additional);
                new_d = zeros(additional * 2, 1);

                predIds = (1:additional)';
                
                rows = [predIds * 2 - 1; predIds * 2];

                cols = [current_nVar + predIds; current_nVar + predIds];

                values = [zeros(additional, 1) - 1; zeros(additional, 1) + 1];

                new_C(sub2ind(size(new_C),rows, cols)) = values;

                new_d(predIds * 2 - 1) = 1;
                new_d(predIds * 2) = 1;
                
                new_predicate_lb = [In.predicate_lb; zeros(additional, 1)];
                new_predicate_ub = [In.predicate_ub; zeros(additional, 1)];

                new_predicate_lb(predIds + current_nVar) = sign_lb;
                new_predicate_ub(predIds + current_nVar) = sign_ub;
                
                new_d = [In.d; new_d];
                new_C = [[In.C zeros(size(In.C, 1), additional)]; new_C];

                if ~isempty(In.Z)
                    new_Z.c = new_V(:, 1);

                    new_Z.V = new_V(:, 2:end);
                end

                S = Star(new_V, new_C, new_d, new_predicate_lb, new_predicate_ub, new_Z);
            end
        end
      
        function S = stepReachApprox_v2(In)
            % @In: star set input
            % @index: index of the neuron performing stepReach
            % @S: star output

            % author: Mykhailo Ivashchenko
            
            [lb, ub] = In.estimateRanges();

            if isempty(lb) || isempty(ub)
                S = [];
            else
                current_nVar = In.nVar;

                l = length(lb);         

                new_V = In.V;
                new_predicate_lb = [In.predicate_lb; zeros(l, 1)];
                new_predicate_ub = [In.predicate_ub; zeros(l, 1)];
                new_Z = In.Z;

                neg_bool = ub < 0;
                pos_bool = lb >= 0;
                if(isa(neg_bool, "dlarray"))
                    neg_bool = extractdata(neg_bool);
                    pos_bool = extractdata(pos_bool);
                end

                neg_ids = find(neg_bool);
                pos_ids = find(pos_bool);

                ids = (1 : length(ub))';
                others = setdiff(ids, neg_ids);
                others = setdiff(others, pos_ids);

                new_V = [new_V eye(l, l)];

                %%% new_V([neg_ids; pos_ids], 2:end) = zeros(length([neg_ids; pos_ids]), size(new_V, 2) - 1);

                new_V(:, 2 : (current_nVar + 1)) = zeros(l, current_nVar);
                new_V(neg_ids, 1) = -1;
                new_V(pos_ids, 1) = 1;
                new_V(others, 1) = Sign.evaluate(new_V(others, 1));

                new_C = zeros((2 * l), current_nVar + l);
                new_d = zeros(l * 2, 1);

                rows_neg1 = neg_ids * 2 - 1;
                rows_neg2 = neg_ids * 2;
                cols = current_nVar + neg_ids;

                new_C(sub2ind(size(new_C),rows_neg1, cols)) = -1;
                new_C(sub2ind(size(new_C),rows_neg2, cols)) = 1;

                rows_pos1 = pos_ids * 2 - 1;
                rows_pos2 = pos_ids * 2;
                cols = current_nVar + pos_ids;

                new_C(sub2ind(size(new_C),rows_pos1, cols)) = -1;
                new_C(sub2ind(size(new_C),rows_pos2, cols)) = 1;

                rows_others1 = others * 2 - 1;
                rows_others2 = others * 2;
                cols = current_nVar + others;

                new_C(sub2ind(size(new_C),rows_others1, cols)) = -1;
                new_C(sub2ind(size(new_C),rows_others2, cols)) = 1;

                new_d(neg_ids * 2 - 1) = 1;
                new_d(neg_ids * 2) = -1;
                new_d(pos_ids * 2 - 1) = -1;
                new_d(pos_ids * 2) = 1;
                new_d(others * 2 - 1) = 1;
                new_d(others * 2) = 1;

                new_predicate_lb(neg_ids + current_nVar) = -1;
                new_predicate_ub(neg_ids + current_nVar) = -1;
                new_predicate_lb(pos_ids + current_nVar) = 1;
                new_predicate_ub(pos_ids + current_nVar) = 1;
                new_predicate_lb(others + current_nVar) = -1;
                new_predicate_ub(others + current_nVar) = 1;
                
                new_d = [In.d; new_d];
                new_C = [[In.C zeros(size(In.C, 1), l)]; new_C];
                
                if ~isempty(In.Z)
                    new_Z.c = new_V(:, 1);

                    new_Z.V = new_V(:, 2:end);
                end

                S = Star(new_V, new_C, new_d, new_predicate_lb, new_predicate_ub, new_Z);
            end
        end
      
        function S = stepReachApprox_v3(In)
            % @In: star set input
            % @index: index of the neuron performing stepReach
            % @S: star output

            % author: Mykhailo Ivashchenko
            
            [lb, ub] = In.estimateRanges();
            
            if isempty(lb) || isempty(ub)
                S = [];
            else
                current_nVar = In.nVar;

                l = length(lb);         

                new_V = In.V;
                new_predicate_lb = [In.predicate_lb; zeros(l, 1)];
                new_predicate_ub = [In.predicate_ub; zeros(l, 1)];
                new_Z = In.Z;

                neg_bool = ub < 0;
                pos_bool = lb >= 0;
                if(isa(neg_bool, "dlarray"))
                    neg_bool = extractdata(neg_bool);
                    pos_bool = extractdata(pos_bool);
                end

                neg_ids = find(neg_bool);
                pos_ids = find(pos_bool);

                ids = (1 : length(ub))';
                others = setdiff(ids, neg_ids);
                others = setdiff(others, pos_ids);

                new_V = [new_V eye(l, l)];

                new_V(:, 2 : (current_nVar + 1)) = zeros(l, current_nVar);
                new_V(neg_ids, 1) = -1;
                new_V(pos_ids, 1) = 1;
                new_V(others, 1) = Sign.evaluate(new_V(others, 1));

                new_C = zeros((2 * l), current_nVar + l);
                new_d = zeros(l * 2, 1);

                for i=1:length(neg_ids)
                    new_C(neg_ids(i) * 2 - 1, current_nVar + neg_ids(i)) = -1;
                    new_C(neg_ids(i) * 2, current_nVar + neg_ids(i)) = 1;
                end
                %new_C(neg_ids * 2 - 1, size(In.C, 2) + neg_ids) = -1;
                %new_C(neg_ids * 2, size(In.C, 2) + neg_ids) = 1;
                new_d(neg_ids * 2 - 1) = 1;
                new_d(neg_ids * 2) = -1;

                new_predicate_lb(neg_ids + current_nVar) = -1;
                new_predicate_ub(neg_ids + current_nVar) = -1;

                for i=1:length(pos_ids)
                    new_C(pos_ids(i) * 2 - 1, current_nVar + pos_ids(i)) = -1;
                    new_C(pos_ids(i) * 2, current_nVar + pos_ids(i)) = 1;
                end
                %new_C(pos_ids * 2 - 1, size(In.C, 2) + pos_ids) = 1;
                %new_C(pos_ids * 2, size(In.C, 2) + pos_ids) = -1;
                new_d(pos_ids * 2 - 1) = -1;
                new_d(pos_ids * 2) = 1;

                new_predicate_lb(pos_ids + current_nVar) = 1;
                new_predicate_ub(pos_ids + current_nVar) = 1;

                for i=1:length(others)
                    new_C(others(i) * 2 - 1, current_nVar + others(i)) = -1;
                    new_C(others(i) * 2, current_nVar + others(i)) = 1;
                end
                %new_C(others * 2 - 1, size(In.C, 2) + others) = 1;
                %new_C(others * 2, size(In.C, 2) + others) = -1;
                new_d(others * 2 - 1) = 1;
                new_d(others * 2) = 1;

                new_predicate_lb(others + current_nVar) = -1;
                new_predicate_ub(others + current_nVar) = 1;

                new_d = [In.d; new_d];
                new_C = [[In.C zeros(size(In.C, 1), l)]; new_C];
                
                if ~isempty(In.Z)
                    new_Z.c = new_V(:, 1);

                    new_Z.V = new_V(:, 2:end);
                end

                S = Star(new_V, new_C, new_d, new_predicate_lb, new_predicate_ub, new_Z);
            end
        end   
        
        % over-approximate reachability analysis using Star
        function S = reach_star_approx(I, mode)
            % @I: star input set
            % @S: star output set

            % author: Mykhailo Ivashchenko
            
            
            if isempty(I)
                S = [];
            else                
                if isa(I, 'ImageStar')
                    Converted_star = I.toStar();
                    In = Star(Converted_star.V, Converted_star.C, Converted_star.d, Converted_star.predicate_lb, Converted_star.predicate_ub, Converted_star.Z);
                else
                    In = Star(I.V, I.C, I.d, I.predicate_lb, I.predicate_ub, I.Z);
                end
                S = Sign.stepReachApprox_v1(In, mode);
                                    
                if isa(I, 'ImageStar')
                    S = S.toImageStar(I.height, I.width, I.numChannel);
                end
            end
                        
        end    
    end
    
    methods(Static) % main reach method
        
        % main function for reachability analysis
        function R = reach(varargin)
            % @I: an array of star input sets
            % @method: 'exact-star' or 'approx-star'
            % @option: = 'parallel' use parallel option
            %          = '' do use parallel option
            
            % author: Mykhailo Ivashchenko
            
            switch nargin
                
                case 7
                    I = varargin{1};
                    method = varargin{2};
                    option = varargin{3};
                    relaxFactor = varargin{4}; % used for aprox-star only
                    dis_opt = varargin{5}; % display option
                    lp_solver = varargin{6}; 
                    mode = varargin{7};
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
                R = Sign.reach_star_exact_multipleInputs(I, option, dis_opt, lp_solver, mode);
            elseif strcmp(method, 'approx-star')  % over-approximate analysis using star
                R = Sign.reach_star_approx(I, mode);
            % TODO: add "not implemented" for other methods
            end
                            
        end
        
        
    end
    
    
end