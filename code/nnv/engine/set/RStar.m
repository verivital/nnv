classdef RStar
    % Relaxed Star class
    % Sung Woo Choi: 12/18/2020
    
    properties
        V = []; % basic matrix that contains c and X
        C = []; % constraint matrix
        d = []; % constraint vector
        c = []; % center vector from V
        X = []; % basic matrix from V
        
        predicate_lb = []; % lower bound vector of predicate variable
        predicate_ub = []; % upper bound vector of predicate variable
        
        lower_a = {[]}; % a set of matrix for lower constraint for bound (a[<=]) ((1 a[1] a[2] ... a[n])'
        upper_a = {[]}; % a set of matrix for upper constraint for bound (a[>=])
        lb = {[]}; % a set of matrix for lower bound
        ub = {[]}; % a set of matrix for upper bound

        iter = inf; % number of iterations for back substitution
        dim = 0; % dimension of current relaxed star set
    end
    
    methods
        
        % constructor
        function obj = RStar(varargin)
            % @V: bassic matrix
            % @C: constraint matrix
            % @d: constraint vector
            % @lower_a: a set of matrix for lower polyhedral constraint
            % @upper_a: a set of matrix for upper polyhedral constraint
                
            switch nargin
                
                case 10
                    V = varargin{1};
                    C = varargin{2};
                    d = varargin{3};
                    pred_lb = varargin{4};
                    pred_ub = varargin{5};
                    lower_a = varargin{6};
                    upper_a = varargin{7};
                    lb = varargin{8};
                    ub = varargin{9};
                    iter = varargin{10};
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % add code for checking properties
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    [nV, mV] = size(V);
                    [nC, mC] = size(C);
                    [nd, md] = size(d);
                    [n1, m1] = size(pred_lb);
                    [n2, m2] = size(pred_ub);
                    
                    if mV ~= mC + 1
                        error('Inconsistency between basic matrix and constraint matrix');
                    end

                    if nC ~= nd
                        error('Inconsistency between constraint matrix and constraint vector');
                    end

                    if md ~= 1
                        error('constraint vector should have one column');
                    end
                    
                    if m1 ~=1 || m2 ~=1 
                        error('predicate lower- or upper-bounds vector should have one column');
                    end
                    
                    if n1 ~= n2 || n1 ~= mC
                        error('Inconsistency between number of predicate variables and predicate lower- or upper-bounds vector');
                    end
                    
                    if ~iscell(lower_a) || ~iscell(upper_a) || ~iscell(lb) || ~iscell(ub)
                       error('polyhedral constraints and bound must be cell array');
                    end
                    obj.V = V;
                    obj.C = C;
                    obj.d = d;
                    obj.predicate_lb = pred_lb;
                    obj.predicate_ub = pred_ub;
                    
                    obj.lower_a = lower_a;
                    obj.upper_a = upper_a;
                    obj.lb = lb;
                    obj.ub = ub;
                    
                    obj.iter = iter;
                    len = length(lower_a);
                    obj.dim = size(lower_a{len}, 1);
                    
                case 9
                    V = varargin{1};
                    C = varargin{2};
                    d = varargin{3};
                    pred_lb = varargin{4};
                    pred_ub = varargin{5}
                    lower_a = varargin{6};
                    upper_a = varargin{7};
                    lb = varargin{8};
                    ub = varargin{9};
                    iter = inf;
                    
                    obj = RStar(V, C, d, pred_lb, pred_ub, lower_a, upper_a, lb, ub, iter);
                    
                case 3
                    % construct rstar from lower and upper bound vector
                    lb = varargin{1};
                    ub = varargin{2};
                    iter = varargin{3};

                    S = Star(lb, ub);
                    obj.V = S.V;
                    obj.C = S.C;    %zeros(1, S.nVar); % initiate an obvious constraint
                    obj.d = S.d;    %zeros(1, 1);
                    obj.predicate_lb = -ones(S.nVar, 1);
                    obj.predicate_ub = ones(S.nVar, 1);
                    
                    dim = size(lb, 1);                   
                    obj.lower_a{1} = [zeros(dim, 1) -eye(dim)]; 
                    obj.upper_a{1} = [zeros(dim, 1) eye(dim)];
                    obj.lb{1} = lb;
                    obj.ub{1} = ub;
                    
                    obj.iter = iter;
                    obj.dim = dim;
                case 2
                    I = varargin{1};
                    iter = varargin{2};
                    if isa(I, 'Polyhedron')
                        dim = I.Dim;
                        c = zeros(dim,1);
                        Ve = eye(dim);
                        V = [c Ve];
                        if ~isempty(I.Ae)
                            C = [I.A;I.Ae;-I.Ae];
                            d = [I.b;I.be;-I.be];
                        else
                            C = I.A;
                            d = I.b;
                        end
                        m = size(C,2);
                        pred_lb = -ones(m, 1);
                        pred_ub = ones(m, 1);
                        I.outerApprox;
                        l = I.Internal.lb;
                        u = I.Internal.ub;
                    elseif isa(I, 'Star')
                        dim = I.dim;
                        V = I.V;
                        C = I.C;
                        d = I.d;
                        pred_lb = I.predicate_lb;
                        pred_ub = I.predicate_ub;
                        [l, u] = I.getRanges();
                    elseif isa(I, 'Zono')
                        dim = I.dim;
                        [l, u] = I.getRanges();
                        V = [I.c I.V];
                        C = [eye(dim); -eye(dim)];
                        d = [u; -l];
                        m = size(C,2);
                        pred_lb = -ones(m, 1);
                        pred_ub = ones(m, 1);
                    else
                        error('Unkown input set');
                    end
                    
                    if iter <= 0
                        error('Iteration must be greater than zero');
                    end
                    
                    
                    lower_a{1} = [zeros(dim, 1) eye(dim)];
                    upper_a{1} = [zeros(dim, 1) eye(dim)];
                    lb{1} = l;
                    ub{1} = u;
                    obj = RStar(V, C, d, pred_lb, pred_ub, lower_a, upper_a, lb, ub, iter);
                    
                case 1
                    I = varargin{1};
                    n = length(I);
                    if n > 1
                        for i = 1:n
                            obj(i) = RStar(I(i), inf);
                        end
                    else
                        obj = RStar(I, inf);
                    end
                    
                case 0
                    obj.V = [];
                    obj.C = [];
                    obj.d = [];
                    obj.predicate_lb = [];
                    obj.predicate_ub = [];
                    obj.lower_a = {[]};
                    obj.upper_a = {[]};
                    obj.lb = {[]};
                    obj.ub = {[]};
                    obj.iter = inf;
                    obj.dim = 0;
                    
                otherwise
                    error('Invalid number of input arguments (should be 0, 1, 2, 7, 8)');
            end
        end
        
        % intersection with a half space: H(x) := Hx <= g
        function R = intersectHalfSpace(obj, H, g)
            % @H: HalfSpace matrix
            % @g: HalfSpace vector
            % @R: new rstar set with more constraints
            
            % author: Sung Woo Choi
            % date: 12/23/2020
            % reference: star set
            
            [nH, mH] = size(H);
            [ng, mg] = size(g);
            
            if mg ~= 1
                error('Halfspace vector should have one column');
            end
            if nH ~= ng
                error('Inconsistent dimension between Halfspace matrix and halfspace vector');
            end
            if mH ~= obj.dim
                error('Inconsistent dimension between halfspace and star set');
            end
            
            m = size(obj.V, 2);
            
            C1 = H*obj.V(:, 2:m);
            d1 = g - H*obj.V(:,1);
            
            new_C = vertcat(obj.C, C1);
            new_d = vertcat(obj.d, d1);
            
            R = RStar(obj.V, new_C, new_d, obj.predicate_lb, obj.predicate_ub, obj.lower_a, obj.upper_a, obj.lb, obj.ub, obj.iter);
            if R.isEmptySet
                R = [];
            end 
        end
        
        % affine abstract mapping of RStar set
        function R = affineMap(varargin)
            switch nargin
                case 3
                    obj = varargin{1};
                    W = varargin{2};
                    b = varargin{3};
                case 2
                    obj = varargin{1};
                    W = varargin{2};
                    b = zeros(size(W,1),1);
            end
            
            [nW, mW] = size(W);
            [nb, mb] = size(b);
            
            if mW ~= obj.dim
                error('Inconsistency between the affine mapping matrix and dimension of the RStar set');
            end
            
            if mb > 1
                error('bias vector must be one column');
            end
            
            if nW ~= nb && nb~=0
                error('Inconsistency between the affine mapping matrix and the bias vector');
            end

            if nb == 0
                b = zeros(nW, 1);
            end

            % affine mapping of basic matrix
            V = W * obj.V;
            if mb ~= 0
                V(:, 1) = V(:, 1) + b;
            end

            % new lower and upper polyhedral contraints
            lower_a = obj.lower_a;
            upper_a = obj.upper_a;
            len = length(lower_a);
            
            lower_a{len+1} = [b W];
            upper_a{len+1} = [b W];
            
            % new lower and uppper bounds
            lb = obj.lb;
            ub = obj.ub;
            lb{len+1} = obj.lb_backSub(lower_a, upper_a);
            ub{len+1} = obj.ub_backSub(lower_a, upper_a);

            R = RStar(V, obj.C, obj.d, obj.predicate_lb, obj.predicate_ub, lower_a, upper_a, lb, ub, obj.iter); 
        end

        % lower bound back-substitution
        function lb = lb_backSub(obj, lower_a, upper_a)
            maxIter = obj.iter;
            len = length(upper_a);
            [nL, mL] = size(upper_a{len});
            alpha = upper_a{len}(:,2:end);
            lower_v = zeros(nL, 1);
            upper_v = upper_a{len}(:,1);
            
            % b[s+1] = v' + sum( max(0,w[j]')*lower_a[j] + min(w[j]',0)*upper_a[j}] ) for j is element of k and for k < i
            % iteration until lb' = b[s'] = v''
            len = len - 1;
            iter = 0;
            while (len > 1 && iter < maxIter)
                [nL, mL] = size(upper_a{len});
                dim = nL;
                
                max_a = max(0, alpha);
                min_a = min(alpha, 0);

                lower_v = max_a * lower_a{len}(:,1) + lower_v;
                upper_v = min_a * upper_a{len}(:,1) + upper_v;
                
                alpha = max_a * lower_a{len}(:,2:end) + ...
                        min_a * upper_a{len}(:,2:end);

                len = len - 1;
                iter = iter + 1;
            end
            
            max_a = max(0, alpha);
            min_a = min(alpha, 0);
            
            [lb1,ub1] = getRanges_L(obj,len);
            lb = max_a * lb1 + lower_v + ...
                 min_a * ub1 + upper_v;
        end
        
        % upper bound back-substituion
        function ub = ub_backSub(obj, lower_a, upper_a)
            maxIter = obj.iter;
            len = length(upper_a);
            [nL, mL] = size(upper_a{len});
            alpha = upper_a{len}(:,2:end);
            lower_v = zeros(nL, 1);
            upper_v = upper_a{len}(:,1);
            
            % c[t+1] = v' + sum( max(0,w[j]')*upper_a[j] + min(w[j]',0)*lower_a[j}] )  for j is element of k and for k < i
            % iteration until ub' = c[t'] = v''
            len = len - 1;
            iter = 0;
            while (len > 1 && iter < maxIter)
                dim = size(lower_a{len}, 1);
                
                max_a = max(0, alpha);
                min_a = min(alpha, 0);
                
                lower_v = min_a * lower_a{len}(:,1) + lower_v;
                upper_v = max_a * upper_a{len}(:,1) + upper_v;
                
                alpha = min_a * lower_a{len}(:,2:end) + ...
                        max_a * upper_a{len}(:,2:end);
                    
                len = len - 1;
                iter = iter + 1;
            end
           
            max_a = max(0, alpha);
            min_a = min(alpha, 0);
            
            [lb1,ub1] = getRanges_L(obj,len);
            ub = min_a * lb1 + lower_v + ...
                 max_a * ub1 + upper_v;
        end
        
        % check is empty set
        function bool = isEmptySet(obj)
            
            % author: Sung Woo Choi
            % date: 12/22/2020
            % reference: star set
            
            options = optimoptions(@linprog, 'Display','none'); 
            options.OptimalityTolerance = 1e-10; % set tolerance
            nVar = size(obj.C, 2);
            f = zeros(1, nVar);
            [~, ~, exitflag, ~] = linprog(f, obj.C, obj.d, [], [],  obj.predicate_lb, obj.predicate_ub, options);
            
            if exitflag == 1
                bool = 0;
            elseif exitflag == -2
                bool = 1;
            else
                error('Error, exitflag = %d', exitflag);
            end            
        end
        
        % get exact lower bound and upper bound vector of the state variables
        function [lb, ub] = getExactRanges(obj)
            
            % author: Sung Woo Choi
            % date: 12/22/2020
            % reference: star set
            
            if ~obj.isEmptySet            
                n = obj.dim;
                lb = zeros(n,1);
                ub = zeros(n,1);
                for i=1:n
                    % fprintf('\nGet range at index %d', i);
                    [lb(i), ub(i)] = obj.getExactRange(i);
                end
            else
                lb = [];
                ub = [];
            end

        end
        
        % get exact lower and upper bounds at specific position using glpk
        function [xmin, xmax] = getExactRange(obj, index)
            % @index: position of the state
            % range: min and max values of x[index]
            
            % author: Sung Woo Choi
            % date: 12/22/2020
            % reference: star set
            
            if index < 1 || index > obj.dim
                error('Invalid index');
            end
            nVar = size(obj.C, 2);
            f = obj.V(index, 2:nVar + 1);
            if all(f(:)==0)
                xmin = obj.V(index,1);
                xmax = obj.V(index,1);
            else
                % **** linprog is much faster than glpk
                options = optimoptions(@linprog, 'Display','none'); 
                options.OptimalityTolerance = 1e-10; % set tolerance
                [~, fval, exitflag, ~] = linprog(f, obj.C, obj.d, [], [], obj.predicate_lb, obj.predicate_ub, options);             
%                 [~, fval, exitflag, ~] = glpk(f, obj.C, obj.d, obj.predicate_lb, obj.predicate_ub);
                if exitflag > 0
                    xmin = fval + obj.V(index, 1);
                else
                    error('Cannot find an optimal solution, exitflag = %d', exitflag);
                end          
          
                [~, fval, exitflag, ~] = linprog(-f, obj.C, obj.d, [], [], obj.predicate_lb, obj.predicate_ub, options);   
%                 [~, fval, exitflag, ~] = glpk(-f, obj.C, obj.d, obj.predicate_lb, obj.predicate_ub);
                if exitflag > 0
                    xmax = -fval + obj.V(index, 1);
                else
                    error('Cannot find an optimal solution');
                end

            end
                        
        end
        
        % get the lower and upper bounds of a current layer at specific
        % position
        function [lb,ub] = getRange(obj, i)
            if i > obj.dim
                error('i should not exceed dimnesion');
            end
            
            len = length(obj.lb);
            lb = obj.lb{len}(i);
            ub = obj.ub{len}(i);
        end
        
        % get lower and upper bounds of a current layer
        function [lb,ub] = getRanges(obj)
            len = length(obj.lb);
            lb = obj.lb{len};
            ub = obj.ub{len};
        end
        
        % get lower and upper bounds of a specific layer
        function [lb,ub] = getRanges_L(obj, len)
            numL = length(obj.lower_a);
            if len > numL
                error('range request should be layers within iteration');
            end
            lb = obj.lb{len};
            ub = obj.ub{len};
        end
        
        % convert to Polyhedron
        function P = toPolyhedron(obj)
            if ~isempty(obj.predicate_lb) && ~isempty(obj.predicate_ub)
                [l, u] = obj.getRanges;
                m = size(obj.C, 2);
                C1 = [eye(m); -eye(m)];
                d1 = [obj.predicate_ub; -obj.predicate_lb];
                
                C = [obj.C; C1];
                d = [obj.d; d1];
            else
                C = obj.C;
                d = obj.d;
            end
            P1 = Polyhedron('A', C, 'b', d);
            X = obj.X;
            c = obj.c;
            P = X*P1 + c;
        end
        
        % convert to Zonotope
        function Z = toZono(obj)
            Z = Zono(obj.c, obj.X);
        end
        
        % convert to Star
        function S = toStar(obj)
            S = Star(obj.V, obj.C, obj.d, obj.predicate_lb, obj.predicate_ub);
        end
        
        % get basic matrix
        function X = get.X(obj)
            X = obj.V(:, 2:end);
        end
        
        % get center vector
        function c = get.c(obj)
            c = obj.V(:,1);
        end
        
        % plot RStar set
        function plot (varargin)                   
            switch nargin
                case 1
                    obj = varargin{1};
                    color = 'red';
                case 2
                    obj = varargin{1};
                    color = varargin{2};
            end
          
            n = length(obj);
            if ~strcmp(color, 'rand')
                c_rand = color;
            end
            
            hold on;
            for i=1:n
                I = obj(i);
                P = I.toPolyhedron;
                if strcmp(color, 'rand')
                    c_rand = rand(1,3);
                end
            
                plot(P, 'color', c_rand);
            end
            hold off
        end
        
        % intersection with other star set (half space)
        function S = Intersect(obj1, obj2)
            C1 = obj2.C * obj1.X;
            d1 = obj2.d - obj2.C * obj1.c;

            new_C = [obj1.C; C1];
            new_d = [obj1.d; d1];
            S = Star(obj1.V, new_C, new_d);     
            if isEmptySet(S)
                S = [];
            end
        end
        
        % check if a index is larger than other
        function bool = is_p1_larger_than_p2(obj, p1_id, p2_id)
            % @p1_id: index of point 1
            % @p2_id: index of point 2
            % @bool = 1 if there exists the case that p1 >= p2
            %       = 0 if there is no case that p1 >= p2
            
            % author: Sung Woo Choi
            % date: 06/23/2021
            % reference: star set
          
            if p1_id < 1 || p1_id > obj.dim
                error('Invalid index for point 1');
            end
            
            if p2_id < 1 || p2_id > obj.dim
                error('Invalid index for point 2');
            end
            
            d1 = obj.c(p1_id) - obj.c(p2_id);
            C1 = obj.X(p2_id, :) - obj.X(p1_id, :);
            S = Star(obj.V, [obj.C; C1], [obj.d;d1], obj.predicate_lb, obj.predicate_ub);
            if S.isEmptySet 
                bool = 0;
            else
                bool = 1;
            end
        end

    end
end