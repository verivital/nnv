classdef Star
    %Star set class 
    %   Star set defined by x = c + a[1]*v[1] + a[2]*v[2] + ... + a[n]*v[n]
    %                         = V * b, V = [c v[1] v[2] ... v[n]]
    %                                  b = [1 a[1] a[2] ... a[n]]^T                                   
    %                       where C*a <= d, constraints on a[i]
    %   Dung Tran: 10/22/2018
    %   updates: Diego Manzanas (01/20/2023)
    %            - lpsolver function substitution
    %            - organized fuctions, clean code and add cast operations
    %                for consistency
    
    properties
        V = []; % basic matrix
        C = []; % constraint matrix
        d = []; % constraint vector
        dim = 0; % dimension of star set
        nVar = 0; % number of variable in the constraints
        predicate_lb = []; % lower bound vector of predicate variable
        predicate_ub = []; % upper bound vector of predicate variable
        state_lb = []; % lower bound of state variables
        state_ub = []; % upper bound of state variables
        Z = []; % an outer zonotope covering this star, used for reachability of logsig and tansig networks
    end
    
    methods % constructor and main methods
        
        % constructor
        function obj = Star(varargin)
            % @V: bassic matrix
            % @C: constraint matrix
            % @d: constraint vector
            
            switch nargin
                
                case 7
                    
                    V = varargin{1};
                    C = varargin{2};
                    d = varargin{3};
                    pred_lb = varargin{4};
                    pred_ub = varargin{5};
                    state_lb = varargin{6};
                    state_ub = varargin{7};
                    
                    [nV, mV] = size(V);
                    [nC, mC] = size(C);
                    [nd, md] = size(d);
                    [n1, m1] = size(pred_lb);
                    [n2, m2] = size(pred_ub);
                    [n3, m3] = size(state_lb);
                    [n4, m4] = size(state_ub);
                    
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
                    
                    if n3 ~= nV || n4 ~= nV
                        error('Inconsistent dimension between lower bound and upper bound vector of state variables and matrix V');
                    end
                    
                    if m3 ~= 1 || m4 ~= 1
                        error('Invalid lower bound or upper bound vector of state variables');
                    end
                        
                    obj.V = V;
                    obj.C = C;
                    obj.d = d;
                    obj.dim = nV;
                    obj.nVar = mC;
                    obj.predicate_lb = pred_lb;
                    obj.predicate_ub = pred_ub;
                    obj.state_lb = state_lb;
                    obj.state_ub = state_ub;
                    
                case 6
                    
                    V = varargin{1};
                    C = varargin{2};
                    d = varargin{3};
                    pred_lb = varargin{4};
                    pred_ub = varargin{5};
                    outer_zono = varargin{6};
                    
                    [nV, mV] = size(V);
                    [nC, mC] = size(C);
                    [nd, md] = size(d);
                    
                    if ~isempty(outer_zono) && ~isa(outer_zono, 'Zono')
                        error("Outer zonotope is not a Zono object");
                    end
                    
                    if ~isempty(outer_zono)
                        [nZ, ~] = size(outer_zono.V);
                        if nZ ~= nV
                            error('Inconsistent dimension between outer zonotope and star set');
                        end
                    end

                    if mV ~= mC + 1
                        error('Inconsistency between basic matrix and constraint matrix');
                    end
                    
                    if nC ~= nd
                        error('Inconsistency between constraint matrix and constraint vector');
                    end

                    if md ~= 1
                        error('constraint vector should have one column');
                    end
                    
                    if ~isempty(pred_lb) && ~isempty(pred_ub)
                        [n1, m1] = size(pred_lb);
                        [n2, m2] = size(pred_ub);
                        
                        if m1 ~=1 || m2 ~=1 
                            error('predicate lower- or upper-bounds vector should have one column');
                        end

                        if n1 ~= n2 || n1 ~= mC
                            error('Inconsistency between number of predicate variables and predicate lower- or upper-bounds vector');
                        end
                    end
                    
                    obj.V = V;
                    obj.C = C;
                    obj.d = d;
                    
                    obj.dim = nV;
                    obj.nVar = mC;
                    obj.predicate_lb = pred_lb;
                    obj.predicate_ub = pred_ub;
                    obj.Z = outer_zono;
                
                case 5
                    
                    V = varargin{1};
                    C = varargin{2};
                    d = varargin{3};
                    pred_lb = varargin{4};
                    pred_ub = varargin{5};
                    
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
                    
                    if (m1 ~= 0 && m2~=0) && (m1 ~=1 || m2 ~=1) 
                        error('predicate lower- or upper-bounds vector should have one column');
                    end
                    
                    if (n1 ~=0 && n2 ~= 0) && (n1 ~= n2 || n1 ~= mC)
                        error('Inconsistency between number of predicate variables and predicate lower- or upper-bounds vector');
                    end

                    obj.V = V;
                    obj.C = C;
                    obj.d = d;
                    obj.dim = nV;
                    obj.nVar = mC;
                    obj.predicate_lb = pred_lb;
                    obj.predicate_ub = pred_ub;
                
                case 3

                    V = varargin{1};
                    C = varargin{2};
                    d = varargin{3};
                    [nV, mV] = size(V);
                    [nC, mC] = size(C);
                    [nd, md] = size(d);
                    
                    if mV ~= mC + 1
                        error('Inconsistency between basic matrix and constraint matrix');
                    end

                    if nC ~= nd
                        error('Inconsistency between constraint matrix and constraint vector');
                    end

                    if md ~= 1
                        error('constraint vector should have one column');
                    end

                    obj.V = V;
                    obj.C = C;
                    obj.d = d;
                    obj.dim = nV;
                    obj.nVar = mC;
                                                            
                case 2
                    
                    % construct star from lower bound and upper bound vector
                    lb = varargin{1};
                    ub = varargin{2};
                    
                    B = Box(lb,ub);
                    S = B.toStar;
                    obj.V = S.V;
                    obj.C = cast(zeros(1, S.nVar), 'like', lb); % initiate an obvious constraint
                    obj.d = cast(zeros(1, 1), 'like', lb);
                    obj.dim = S.dim;
                    obj.nVar = S.nVar;
                    obj.state_lb = lb;
                    obj.state_ub = ub;
                    obj.predicate_lb = cast(-ones(S.nVar, 1), 'like', lb);
                    obj.predicate_ub = cast(ones(S.nVar, 1), 'like', lb);
                    obj.Z = B.toZono;
                 
                case 1 % accept a polyhedron as an input and transform to a star

                    I = varargin{1};
                    if ~isa(I, 'Polyhedron')
                        error('Input set is not a polyhedron');
                    end
                    
                    c = cast(zeros(I.Dim, 1), 'like', I.A);
                    V1 = cast(eye(I.Dim), 'like', I.A);
                    V = [c V1];
                    if isempty(I.Ae)    
                        obj = Star(V, I.A, I.b);
                    else
                        A1 = [I.Ae; -I.Ae];
                        b1 = [I.be; -I.be];
                        obj = Star(V, [I.A; A1], [I.b; b1]);
                    end
                    [lb, ub] = obj.getRanges;
                    obj.predicate_lb = lb;
                    obj.predicate_ub = ub;
                    B = Box(lb, ub);
                    obj.Z = B.toZono;
                
                case 0

                    % create empty Star (for preallocation an array of star)
                    obj.V = [];
                    obj.C = [];
                    obj.d = [];
                    obj.dim = 0;
                    obj.nVar = 0;
                               
                otherwise
                    error('Invalid number of input arguments (should be 0, 1, 2, 3, 5, 6 or 7)');
            end
            
        end
        
        % sampling a star set (TODO: optimize sampling)
        function V = sample(obj, N)
            % @N: number of points in the samples
            % @x: a set of N sampled points in the star set 
            % author: Dung Tran
            % date: 1/3/2019
            
            if N < 1
                error('Invalid number of samples');
            end
            
            B = obj.getBox;
            if isempty(B)
                V = [];
            else
                lb = B.lb;
                ub = B.ub;
                X = cell(1, obj.dim);
                V1 = [];
                for i=1:obj.dim
                    X{1, i} = (ub(i) - lb(i)).*rand(2*N, 1) + lb(i);
                    V1 = vertcat(V1, X{1, i}');
                end
                V = [];
                for i=1:2*N
                    if obj.contains(V1(:, i))
                        V = [V V1(:, i)];
                    end
                end
                if size(V, 2) > N               
                    V = V(:, 1:N);
                end             
            end
        end
        
        % affine mapping of star set S = Wx + b;
        function S = affineMap(obj, W, b)
            % @W: mapping matrix
            % @b: mapping vector
            % @S: new star set
            
            if size(W, 2) ~= obj.dim
                error('Inconsistency between the affine mapping matrix and dimension of the star set');
            end
            
            if ~isempty(b)
                if size(b, 1) ~= size(W, 1)
                    error('Inconsistency between the mapping vec and mapping matrix');
                end
                if size(b, 2) ~= 1
                    error('Mapping vector should have one column');
                end
                newV = W * obj.V;
                newV(:, 1) = newV(:, 1) + b;
            else
                newV = W * obj.V;
            end
            
            if ~isempty(obj.Z)
                new_Z = obj.Z.affineMap(W, b);
            else
                new_Z = [];
            end
            
            S = Star(newV, obj.C, obj.d, obj.predicate_lb, obj.predicate_ub, new_Z);
                       
        end
        
        % Minkowski Sum
        function S = MinkowskiSum(obj, X)
            % @X: another star with same dimension
            % @S: new star
            
            if ~isa(X, 'Star')
                error('Input matrix is not a Star');
            else
                if X.dim ~= obj.dim
                    error('Input star and current star have different dimensions');
                end
            end
            
            m1 = size(obj.V, 2);
            m2 = size(X.V, 2);
            
            V1 = obj.V(:, 2:m1);
            V2 = X.V(:, 2:m2);
            new_c = obj.V(:, 1) + X.V(:, 1);
            
            % check if two Star has the same constraints
            if size(obj.C, 1) == size(X.C, 1) && size(obj.C, 2) == size(X.C, 2) && norm(obj.C - X.C) + norm(obj.d - X.d) < 0.0001
               V3 = V1 + V2;
               new_V = horzcat(new_c, V3);
               S = Star(new_V, obj.C, obj.d, obj.predicate_lb, obj.predicate_ub); % new Star has the same number of basic vectors
            else
                V3 = horzcat(V1, V2);
                new_c = obj.V(:, 1) + X.V(:, 1);
                new_V = horzcat(new_c, V3);        
                new_C = blkdiag(obj.C, X.C);        
                new_d = vertcat(obj.d, X.d);
                new_pred_lb = [obj.predicate_lb; X.predicate_lb];
                new_pred_ub = [obj.predicate_ub; X.predicate_ub];
                S = Star(new_V, new_C, new_d, new_pred_lb, new_pred_ub); % new Star has more number of basic vectors
            end
                
        end
        
        function S = HadamardProduct(obj, X)
            % @X: another star with same dimension
            % @S: new star
            
            if ~isa(X, 'Star')
                error('Input matrix is not a Star');
            else
                if X.dim ~= obj.dim
                    error('Input star and current star have different dimensions');
                end
            end
            
            m1 = size(obj.V, 2);
            m2 = size(X.V, 2);
            V1 = obj.V(:, 2:m1);
            V2 = X.V(:, 2:m2);
            new_c = obj.V(:, 1) .* X.V(:, 1);
            
            % check if two Star has the same constraints
            if size(obj.C, 1) == size(X.C, 1) && size(obj.C, 2) == size(X.C, 2) && norm(obj.C - X.C) + norm(obj.d - X.d) < 0.0001
               V3 = V1 .* V2;
               new_V = horzcat(new_c, V3);
               S = Star(new_V, obj.C, obj.d, obj.predicate_lb, obj.predicate_ub); % new Star has the same number of basic vectors
            else
                V3 = horzcat(V1, V2);
                % new_c = obj.V(:, 1) + X.V(:, 1);
                new_V = horzcat(new_c, V3);        
                new_C = blkdiag(obj.C, X.C);        
                new_d = vertcat(obj.d, X.d);
                new_pred_lb = [obj.predicate_lb; X.predicate_lb];
                new_pred_ub = [obj.predicate_ub; X.predicate_ub];
                S = Star(new_V, new_C, new_d, new_pred_lb, new_pred_ub); % new Star has more number of basic vectors
            end

        end

        % New Minkowski Sum (used for Recurrent Layer reachability)
        function S = Sum(obj, X)
            % @X: another star with same dimension
            % @S: new star
            
            if ~isa(X, 'Star')
                error('Input matrix is not a Star');
            else
                if X.dim ~= obj.dim
                    error('Input star and current star have different dimensions');
                end
            end
            
            m1 = size(obj.V, 2);
            m2 = size(X.V, 2);
            
            V1 = obj.V(:, 2:m1);
            V2 = X.V(:, 2:m2);
            
            V3 = horzcat(V1, V2);
            new_c = obj.V(:, 1) + X.V(:, 1);
            new_V = horzcat(new_c, V3);        
            new_C = blkdiag(obj.C, X.C);        
            new_d = vertcat(obj.d, X.d);
            
            if ~isempty(obj.predicate_lb) && ~isempty(X.predicate_lb)
                new_predicate_lb = [obj.predicate_lb; X.predicate_lb];
                new_predicate_ub = [obj.predicate_ub; X.predicate_ub];
                S = Star(new_V, new_C, new_d, new_predicate_lb, new_predicate_ub);
            else
                S = Star(new_V, new_C, new_d); % new Star has more number of basic vectors
            end
       
        end
        
        % scalar map of a Star S' = alp * S, 0 <= alp <= alp_max
        function S = scalarMap(obj, alp_max)
            % @a_max: maximum value of a
            % @S: new Star
            
            % note: we always require that alp >= 0
            % author: Dung Tran
            % date: 10/27/2018
            
            % =============================================================
            % S: x = alp*c + V* alph * a, Ca <= d
            % note that:   Ca <= d -> C*alph*a <= alp*a <= alp_max * d
            % let: beta = alp * a, we have
            % S := x = alp * c + V * beta, C * beta <= alp_max * d,
            %                              0 <= alp <= alp_max
            % Let g = [beta; alp]
            % S = Star(new_V, new_C, new_d), where:
            %   new_V = [0 c V], new_C = [0 -1; 0 1; 0 C], new_d = [0; alpha_max; alp_max * d]
            %       
            % S has one more basic vector compared with obj
            % =============================================================
            
            new_c = cast(zeros(obj.dim, 1), 'like', obj.V);
            new_V = [obj.V new_c];
            new_C = blkdiag(obj.C, [-1; 1]);
            new_d = vertcat(alp_max*obj.d, 0, alp_max);
            S = Star(new_V, new_C, new_d);

        end
        
        % intersection with a half space: H(x) := Hx <= g
        function S = intersectHalfSpace(obj, H, g)
            % @H: HalfSpace matrix
            % @g: HalfSpace vector
            % @S : new star set with more constraints
            
            % author: Dung Tran
            % date: 10/27/2018
            
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
            
            S = Star(obj.V, new_C, new_d, obj.predicate_lb, obj.predicate_ub);
            
            if S.isEmptySet
                S = [];
            end
            
        end
        
    end


    methods % get methods (also estimate)

        % find a box bounding a star
        function B = getBox(obj)
            
            if isempty(obj.C) || isempty(obj.d) % star set is just a vector (one point)
                lb = obj.V(:, 1);
                ub = obj.V(:, 1);               
                B = Box(lb, ub);

            else % star set is a set
                
                if ~isempty(obj.state_lb) && ~isempty(obj.state_ub)
                    B = Box(obj.state_lb, obj.state_ub);

                else
                    
                    if isa(obj.V, 'single') || isa(obj.C, 'single') || isa(obj.d, 'single')
                        obj = obj.changeVarsPrecision('double');
                    end
                    
                    lb = zeros(obj.dim, 1);
                    ub = zeros(obj.dim, 1);
                    
                    for i=1:obj.dim
                        f = obj.V(i, 2:obj.nVar + 1);
                        if all(f(:)==0)
                            lb(i) = obj.V(i,1);
                            ub(i) = obj.V(i,1);
                        else
                            [fval, exitflag] = lpsolver(f, obj.C, obj.d, [], [], obj.predicate_lb, obj.predicate_ub);
                            if ismember(exitflag, ["l1","g5"]) 
                                lb(i) = fval + obj.V(i, 1);
                            else
                                lb = [];
                                ub = [];
                                break;
                            end
                            [fval, exitflag] = lpsolver(-f, obj.C, obj.d, [], [], obj.predicate_lb, obj.predicate_ub);
                            if ismember(exitflag, ["l1","g5"]) 
                                ub(i) = -fval + obj.V(i, 1);
                            else
                                lb = [];
                                ub = [];
                                break;
                            end
                        end


                    end
                    
                    if isempty(lb) || isempty(ub)
                        B = [];
                    else
                        B = Box(lb, ub);           
                    end
                      
                end
                
            end         
              
        end
    
        % estimates ranges of state vector quickly
        function B = getBoxFast(obj)
            % these ranges are not the exact range, it is an
            % over-approximation of the exact ranges

            % author: Dung Tran
            % date: 5/29/2019
            
            if isempty(obj.C) || isempty(obj.d) % star set is just a vector (one point)
                lb = obj.V(:, 1);
                ub = obj.V(:, 1);
                B = Box(lb, ub);
            else      
                [pred_lb, pred_ub] = obj.getPredicateBounds;         
                B1 = Box(pred_lb, pred_ub);
                Z = B1.toZono;
                Z = Z.affineMap(obj.V(:,2:obj.nVar + 1), obj.V(:,1));
                B = Z.getBox;
            end
                         
        end
        
        % return possible max indexes
        function max_ids = getMaxIndexes(obj)
            % @max_ids: index of the state that can be a max point
            
            % author: Dung Tran
            % date: 4/1/2020
            
            new_rs  = obj.toImageStar(obj.dim, 1, 1);                    
            max_id = new_rs.get_localMax_index([1 1], [obj.dim 1], 1);
            max_ids = max_id(:, 1);

        end

        % get bounds of predicate variables
        function [pred_lb, pred_ub] = getPredicateBounds(obj)
            % author: Dung Tran
            % date: 5/30/2019
            
            if isempty(obj.predicate_lb) || isempty(obj.predicate_ub)
                P = Polyhedron('A', obj.C, 'b', obj.d);
                P.outerApprox;
                pred_lb = P.Internal.lb;
                pred_ub = P.Internal.ub;
            else
                pred_lb = obj.predicate_lb;
                pred_ub = obj.predicate_ub;
            end  
            
        end
        
        % find range of a state at specific position
        function [xmin, xmax] = getRange(obj, index)
            % @index: position of the state
            % range: min and max values of x[index]
            
            % author: Dung Tran
            % date: 11/16/2018
            
            if index < 1 || index > obj.dim
                error('Invalid index');
            end
            
            f = obj.V(index, 2:obj.nVar + 1);
            if all(f(:)==0)
                xmin = obj.V(index,1);
                xmax = obj.V(index,1);
            else
                [fval, exitflag] = lpsolver(f, obj.C, obj.d, [], [], obj.predicate_lb, obj.predicate_ub);
                if ismember(exitflag, ["l1","g5"])
                    xmin = fval + obj.V(index, 1);
                else
                    error("Cannot find an optimal solution, exitflag = " + string(exitflag));
                end
                [fval, exitflag] = lpsolver(-f, obj.C, obj.d, [], [], obj.predicate_lb, obj.predicate_ub);
                if ismember(exitflag, ["l1","g5"])
                    xmax = -fval + obj.V(index, 1);
                else
                    error('Cannot find an optimal solution');
                end

            end
                        
        end
        
        % get min
        function xmin = getMin(varargin)
            % @index: position of the state
            % xmin: min value of x[index]
            
            % author: Dung Tran
            % date: 6/11/2018
            % update: 7/16/2020: add lp_solver option
            
            switch nargin
                case 2
                    obj = varargin{1};
                    index = varargin{2};
                    lp_solver = 'linprog';
                case 3
                    obj = varargin{1};
                    index = varargin{2};
                    lp_solver = varargin{3};
                otherwise
                    error('Invalid number of input arguments, should be 2 or 3');
            end
            
            if index < 1 || index > obj.dim
                error('Invalid index');
            end 
            
            f = obj.V(index, 2:obj.nVar + 1);
            if all(f(:)==0)
                xmin = obj.V(index,1);
            else               
                [fval, exitflag] = lpsolver(f, obj.C, obj.d, [], [], obj.predicate_lb, obj.predicate_ub, lp_solver);
                if ismember(exitflag, ["l1","g5"])  
                    xmin = fval + obj.V(index, 1);
                else
                    error('Cannot find an optimal solution, exitflag = ' + exitflag);
                end
            end
            
        end
        
        % get mins
        function xmin = getMins(varargin)
            % @map: an array of indexes
            % xmin: min values of x[indexes]
            
            % author: Dung Tran
            % date: 7/13/2018
            % update: 7/16/2020: add display option + lp_solver option
            
            switch nargin
                case 5
                    obj = varargin{1};
                    map = varargin{2};
                    par_option = varargin{3};
                    dis_option = varargin{4};
                    lp_solver  = varargin{5};
                case 4
                    obj = varargin{1};
                    map = varargin{2};
                    par_option = varargin{3};
                    dis_option = varargin{4};
                    lp_solver = 'linprog';
                case 3
                    obj = varargin{1};
                    map = varargin{2};
                    par_option = varargin{3};
                    dis_option = [];
                    lp_solver = 'linprog';
                case 2
                    obj = varargin{1};
                    map = varargin{2}; 
                    par_option = 'single';
                    dis_option = [];
                    lp_solver = 'linprog';
                otherwise
                    error('Invalid number of inputs, should be 1, 2, 3, or 4');
            end
            
            n = length(map);
            xmin = cast(zeros(n, 1), 'like', obj.V);
            if isempty(par_option) || strcmp(par_option, 'single') % get Maxs using single core
                reverseStr = '';
                for i = 1:n
                    xmin(i) = obj.getMin(map(i), lp_solver);
                    if strcmp(dis_option, 'display')
                        msg = sprintf('%d/%d', i, n);
                        fprintf([reverseStr, msg]);
                        reverseStr = repmat(sprintf('\b'), 1, length(msg));
                    end
                end
            elseif strcmp(par_option, 'parallel') % get Maxs using multiple cores 
                f = obj.V(map, 2:obj.nVar + 1);
                V1 = obj.V(map, 1);
                C1 = obj.C;
                d1 = obj.d;
                pred_lb = obj.predicate_lb;
                pred_ub = obj.predicate_ub;
                parfor i=1:n                    
                    if all(f(i,:)==0)
                        xmin(i) = V1(i,1);
                    else
                        [fval, exitflag] = lpsolver(f(i, :), C1, d1, [], [], pred_lb, pred_ub, lp_solver);
                        if ismember(exitflag, ["l1","g5"])
                            xmin(i) = fval + V1(i, 1);
                        else
                            error('Cannot find an optimal solution, exitflag = %d', exitflag);
                        end                                   
                    end
                end
            else
                error('Unknown parallel option');
            end
            
        end
                
        % get max
        function xmax = getMax(varargin)
            % @index: position of the state
            % xmax: max value of x[index]
            
            % author: Dung Tran
            % date: 6/11/2018
            
             switch nargin
                case 2
                    obj = varargin{1};
                    index = varargin{2};
                    lp_solver = 'linprog';
                case 3
                    obj = varargin{1};
                    index = varargin{2};
                    lp_solver = varargin{3};
                otherwise
                    error('Invalid number of input arguments, should be 2 or 3');
            end
            
            if index < 1 || index > obj.dim
                error('Invalid index');
            end
            
            f = obj.V(index, 2:obj.nVar + 1);
            if all(f(:)==0)
                xmax = obj.V(index,1);
            else
                [fval, exitflag] = lpsolver(-f, obj.C, obj.d, [], [], obj.predicate_lb, obj.predicate_ub, lp_solver);
                if ismember(exitflag, ["l1","g5"])
                    xmax = -fval + obj.V(index, 1);
                else
                    error('Cannot find an optimal solution, exitflag = %s', exitflag);
                end      
            end
            
        end
        
        % get maxs
        function xmax = getMaxs(varargin)
            % @map: an array of indexes
            % xmax: max values of x[indexes]
            
            % author: Dung Tran
            % date: 6/11/2018
            
            switch nargin
                case 5
                    obj = varargin{1};
                    map = varargin{2};
                    par_option = varargin{3};
                    dis_option = varargin{4};
                    lp_solver  = varargin{5};
                case 4
                    obj = varargin{1};
                    map = varargin{2};
                    par_option = varargin{3};
                    dis_option = varargin{4};
                    lp_solver = 'linprog';
                case 3
                    obj = varargin{1};
                    map = varargin{2};
                    par_option = varargin{3};
                    dis_option = [];
                    lp_solver = 'linprog';
                case 2
                    obj = varargin{1};
                    map = varargin{2}; 
                    par_option = 'single';
                    dis_option = [];
                    lp_solver = 'linprog';
                otherwise
                    error('Invalid number of inputs, should be 1, 2, 3, or 4');
            end
            
            n = length(map);
            xmax = cast(zeros(n, 1), 'like', obj.V);
                        
            if isempty(par_option) || strcmp(par_option, 'single') % get Maxs using single core
                reverseStr = '';
                for i = 1:n
                    xmax(i) = obj.getMax(map(i), lp_solver);
                    if strcmp(dis_option, 'display')
                        msg = sprintf('%d/%d', i, n);
                        fprintf([reverseStr, msg]);
                        reverseStr = repmat(sprintf('\b'), 1, length(msg));
                    end
                end
            elseif strcmp(par_option, 'parallel') % get Maxs using multiple cores 
                f = obj.V(map, 2:obj.nVar + 1);
                V1 = obj.V(map, 1);
                C1 = obj.C;
                d1 = obj.d;
                pred_lb = obj.predicate_lb;
                pred_ub = obj.predicate_ub;
                parfor i=1:n
                    if all(f(i,:)==0)
                        xmax(i) = V1(i,1);
                    else
                        [fval, exitflag] = lpsolver(-f(i, :), C1, d1, [], [], pred_lb, pred_ub, lp_solver); 
                        if ismember(exitflag, ["l1","g5"])
                            xmax(i) = -fval + V1(i, 1);
                        else
                            error("Cannot find an optimal solution, exitflag = " + string(exitflag));
                        end                                   
                    end
                end
            else
                error('Unknown parallel option');
            end
     
        end

        % get lower bound and upper bound vector of the state variables
        function [lb, ub] = getRanges(obj)
            
            % author: Dung Tran
            % date: 7/19/2019
            
            if ~obj.isEmptySet            
                n = obj.dim;
                lb = zeros(n,1);
                ub = zeros(n,1);
                for i=1:n
                    [lb(i), ub(i)] = obj.getRange(i);
                end
            else
                lb = [];
                ub = [];
            end

        end
        
        % find range of a state at specific position
        function [xmin, xmax] = estimateRange(obj, index)
            % @index: position of the state
            % range: min and max values of x[index]
            
            % author: Dung Tran
            % date: 7/19/2019
            
            if (index < 1) || (index > obj.dim)
                error('Invalid index');
            end
            
            if ~isempty(obj.predicate_lb) && ~isempty(obj.predicate_ub)

                f = obj.V(index, 1:obj.nVar+1);
                xmin = f(1);
                xmax = f(1);

                for i=2:obj.nVar+1
                    if f(i) >= 0
                        xmin = xmin + f(i) * obj.predicate_lb(i-1);
                        xmax = xmax + f(i) * obj.predicate_ub(i-1);
                    else
                        xmin = xmin + f(i) * obj.predicate_ub(i-1);
                        xmax = xmax + f(i) * obj.predicate_lb(i-1);
                    end

                end

            else
               warning('The ranges of predicate variables are unknown to estimate the ranges of the states, we solve LP optimization to get the exact range');
               [xmin, xmax] = obj.getRange(index); 
            end
        end
        
        % estimate range using clip from Stanley Bak
        function [xmin, xmax] = estimateBound(obj, index)
            % @index: position of the state
            % range: min and max values of x[index]
            
            % author: Dung Tran
            % date: 3/27/2020
            
            if (index < 1) || (index > obj.dim)
                error('Invalid index');
            end
            
            f = obj.V(index, 2:obj.nVar+1);
                        
            pos_mat = f; 
            neg_mat = f;
            pos_mat(pos_mat < 0) = 0; 
            neg_mat(neg_mat > 0) = 0;
            
            xmin1 = pos_mat*obj.predicate_lb;
            xmax1 = pos_mat*obj.predicate_ub;
            xmin2 = neg_mat*obj.predicate_ub;
            xmax2 = neg_mat*obj.predicate_lb;
            
            xmin = obj.V(index, 1) + xmin1 + xmin2;
            xmax = obj.V(index, 1) + xmax1 + xmax2;
            
        end
        
        % quickly estimate lower bound and upper bound vector of statevariables
        function [lb, ub] = estimateBounds(obj)
            
            % author: Dung Tran
            % date: 7/19/2019
            % update: 4/2/2020
            
            if ~isempty(obj.Z)
                [lb, ub] = obj.Z.getBounds;
            else
                n = obj.dim;
                lb = zeros(n,1);
                ub = zeros(n,1);

                for i=1:n
                    [lb(i), ub(i)] = obj.estimateRange(i);
                end
            end

        end
        
        % estimate ranges using clip method from Stanley Bak (slower than for-loop method)
        function [lb, ub] = estimateRanges(obj)
            % @lb: lowerbound vector
            % @ub: upper bound vector
            
            % author: Dung Tran
            % date: 3/27/2020
                                   
            pos_mat = obj.V; 
            neg_mat = obj.V;
            pos_mat(pos_mat < 0) = 0; 
            neg_mat(neg_mat > 0) = 0;

            % ensure predicate bounds are not empty (if empty, set to 1 and -1)
            if isempty(obj.predicate_lb)
                obj.predicate_lb = -1*ones(obj.nVar,1); 
            end
            if isempty(obj.predicate_ub)
                obj.predicate_ub = ones(obj.nVar,1);
            end
            
            xmin1 = pos_mat*[0; obj.predicate_lb];
            xmax1 = pos_mat*[0; obj.predicate_ub];
            xmin2 = neg_mat*[0; obj.predicate_ub];
            xmax2 = neg_mat*[0; obj.predicate_lb];
            
            lb = obj.V(:, 1) + xmin1 + xmin2;
            ub = obj.V(:, 1) + xmax1 + xmax2;
            
        end
        
        % estimate quickly max-point candidates
        function max_cands = get_max_point_candidates(obj)
            %@max_cands: an array of indexes of max-point candidates
            
            % author: Dung Tran
            % date: 7/29/2021
            
            [lb, ub] = obj.estimateRanges();
            [a, id] = max(lb);
            b = (ub >= a);
            if sum(b) == 1
                max_cands = id;
            else
                max_cands = find(b);
            end
                   
        end
        
        % find a oriented box bounding a star
        function B = getOrientedBox(obj)
            % author: Dung Tran
            % date: 11/8/2018
            
            [Q, Z, P] = svd(obj.V(:, 2:obj.nVar + 1));
            S = Z * P';
            lb = cast(zeros(obj.dim, 1), 'like', obj.V);
            ub = cast(zeros(obj.dim, 1), 'like', obj.V);

            for i=1:obj.dim
                f = S(i, :);
                [fval, exitflag] = lpsolver(f, obj.C, obj.d, [], [], obj.predicate_lb, obj.predicate_ub);
                
                if ismember(exitflag, ["l1","g5"])
                    lb(i) = fval;
                else
                    error('Cannot find an optimal solution, exitflag = %d', exitflag);
                end
 
                [fval, exitflag] = lpsolver(-f, obj.C, obj.d, [], [], obj.predicate_lb, obj.predicate_ub);
                if ismember(exitflag, ["l1","g5"])
                    ub(i) = -fval;
                else
                    error('Cannot find an optimal solution, exitflag = %d', exitflag);
                end
            end
            
            new_V = [obj.V(:,1) Q];
            new_C = cast(vertcat(eye(obj.dim), -eye(obj.dim)), 'like', obj.V);
            new_d = vertcat(ub, -lb);
            
            B = Star(new_V, new_C, new_d);
            
        end
        
        % find a zonotope bounding a star (an over-approximation of a star using zonotope)
        function Z = getZono(obj)
            
            % author: Dung Tran
            % date: 10/25/2018
            % update: Diego Manzanas Lopez
            %      - get stored Star.Z if available
            
            if isempty(obj.Z)
                B = obj.getBox;
                if ~isempty(B)
                    Z = B.toZono;
                else
                    Z = [];
                end
            else
                Z = obj.Z;
            end
                        
        end
        
    end


    methods % check methods

        % check if a star set contain a point
        function bool = contains(obj, s)
            % @s: a star point
            % @bool: = 1 star set contains s, else no
            
            % author: Dung Tran
            % date: 1/3/2019
            
            
            if size(s,1) ~= obj.dim
                error('Dimension mismatch');
            end
            if size(s,2) ~= 1
                error('Invalid star point');
            end
            
            A = obj.C;
            b = obj.d;
            Ae = obj.V(:, 2:obj.nVar+1);
            be = s - obj.V(:,1);
            
            P = Polyhedron('A', A, 'b', b, 'Ae', Ae, 'be', be, 'lb', obj.predicate_lb, 'ub', obj.predicate_ub);
            
            bool = ~P.isEmptySet;
                     
        end
        
        % check is empty set
        function bool = isEmptySet(obj)
            % author: Dung Tran
            % date: 
            % update: 6/16/2020
            % update: 7/15/2020 The isEmptySet method in Polyhedron object has bug
            
            f = zeros(1, obj.nVar, 'like', obj.V); % objective function
            [~, exitflag] = lpsolver(f, obj.C, obj.d, [], [], obj.predicate_lb, obj.predicate_ub, 'linprog', 'emptySet');
            if ismember(exitflag,["l1", "g2", "g5"])
                bool = 0;
            elseif ismember(exitflag, ["l-2", "l-5", "g3", "g4", "g110"])
                bool = 1;
            else
                error('Error, exitflag = %d', exitflag);
            end
            
        end
        
        % check if a star set is a subset of other
        function bool = isSubSet(obj, S)
            
            if ~isa(S, 'Star')
                error('Input set is not a Star');
            end
            
            if obj.dim ~= S.dim
                error('Two sets have different dimension');
            end

            if isa(obj.V, 'single') || isa(obj.C,'single') || isa(obj.d, 'single')
                S1 = obj.changeVarsPrecision('double');
            else
                S1 = obj;
            end

            if isa(S.V, 'single') || isa(S.C,'single') || isa(S.d, 'single')
                S2 = S.changeVarsPrecision('double');
            else
                S2 = S;
            end
            
            P1 = S1.toPolyhedron();
            P2 = S2.toPolyhedron();
            
            bool = (P1 <= P2);
            
        end
        
        % check if a index is larger than other
        function bool = is_p1_larger_than_p2(obj, p1_id, p2_id)
            % @p1_id: index of point 1
            % @p2_id: index of point 2
            % @bool = 1 if there exists the case that p1 >= p2
            %       = 0 if there is no case that p1 >= p2
            
            % author: Dung Tran
            % date: 7/10/2020
            
            
            if p1_id < 1 || p1_id > obj.dim
                error('Invalid index for point 1');
            end
            
            if p2_id < 1 || p2_id > obj.dim
                error('Invalid index for point 2');
            end
                        
            d1 = obj.V(p1_id, 1) - obj.V(p2_id, 1); 
            C1 = obj.V(p2_id, 2: obj.nVar + 1) - obj.V(p1_id, 2:obj.nVar+1);
            S = Star(obj.V, [obj.C; C1], [obj.d;d1], obj.predicate_lb, obj.predicate_ub);
            if S.isEmptySet 
                bool = 0;
            else
                bool = 1;
            end 
            
        end
    
    end


    methods % transform / conversion methods
        
        % convex hull of two Stars
        function S = convexHull(obj, X)
            % @X: input star
            % @S: an over-approximation of (convex hull) of two Stars
            
            % author: Dung Tran
            % date: 10/27/2018
            
            % =============================================================
            % This method is similar to the method proposed by Prof. Girard
            % for Zonotope. See the following paper:
            % 1) Reachability of Uncertain Linear Systems Using
            % Zonotopes, Antoin Girard, HSCC2005.
            %
            % CH(S1 U S2) := {l*x1 + (1-l)*x2}, 0 <= l <= 1, x1 \in S1, x2
            % \in S2.     := {x2 + l*(x1 - x2)}
            % 
            %             := c2 + a2 * V2 + l(c1 - c2) + l*a1*V1 - l*a2*V2
            % let beta1 = l * a1, beta2 = l*a2 -> C1 * beta1 <= l * d1 <=
            % d1 and C2 * beta2 <= l * d2 <= d2
            %
            % CH(S1 U S2) = c2 + a2 * V2 + l * (c1 - c2) + beta1 * V1 -
            %                                                   beta2 * V2
            %
            % if S2 = L*S1, i.e., S2 is a projection of S1
            % 
            % CH(S1 U S2) := {(l + (1-l)*L) * x1} = {(l(I - L) + L) * x1}
            %             := (l * (I - L) * x1 + L * x1)
            % =============================================================
            
            if ~isa(X, 'Star')
                error('Input set is not a Star');
            end
            
            if X.dim ~= obj.dim
                error('Inconsisten dimension between input set and this star');
            end
            
            c1 = obj.V(:, 1);
            c2 = X.V(:, 1);
            V1 = obj.V(:, 2:obj.nVar + 1);
            V2 = X.V(:, 2:obj.nVar + 1);
            
            new_V = horzcat(c1-c2, V1, -V2);
            new_C = blkdiag(obj.C, X.C);
            new_d = vertcat(obj.d, X.d);
            
            S2 = Star(new_V, new_C, new_d);
            S = X.MinkowskiSum(S2);
                                    
        end
        
        % convexHull of a Star with its linear transformation
        function S = convexHull_with_linearTransform(obj, L)
            % @L: linear transformation matrix
            % @S: new Star
            
            % author: Dung Tran
            % date: 10/27/2018
            
            [nL, mL] = size(L);
            
            if nL ~= mL 
                error('Trasformation matrix should be a square matrix');
            end
            if nL ~= obj.dim
                error('Inconsistent dimension of tranformation matrix and this zonotope');
            end
            
            M = eye(obj.dim) - L; 
            S1 = obj.affineMap(L, []);
            S2 = obj.affineMap(M, []);
            
            new_V = [S1.V(:,1) S1.V(:,2:S1.nVar+1) S2.V(:,1) S2.V(:,2:S2.nVar+1)];
            new_C = blkdiag(S1.C, [-1; 1], S2.C);
            new_d = vertcat(S1.d, [0; 1], S2.d);
            
            S = Star(new_V, new_C, new_d);
                        
        end
                
        % convert to polyhedron
        function P = toPolyhedron(obj)
            
            b = obj.V(:, 1);        
            W = obj.V(:, 2:obj.nVar + 1);
            
            if ~isempty(obj.predicate_ub)
                C1 = cast([eye(obj.nVar); -eye(obj.nVar)], 'like', obj.V);
                d1 = [obj.predicate_ub; -obj.predicate_lb];
                Pa = Polyhedron('A', [obj.C;C1], 'b', [obj.d;d1]);
                P = W*Pa + b;
            else
                Pa = Polyhedron('A', [obj.C], 'b', [obj.d]);
                P = W*Pa + b;
            end
        end
        
        % convert to ImageStar set
        function IS = toImageStar(obj, height, width, numChannel)
            % @height: height of ImageStar
            % @width: width of ImageStar
            % @numChannel: number of channels in ImageStar
            
            % author: Dung Tran
            % date: 7/19/2019
            
            if obj.dim ~= height*width*numChannel
                error('Inconsistent dimension in the ImageStar and the original Star set');
            end
            
            if ~isempty(obj.state_lb) && ~isempty(obj.state_ub)
                lb = reshape(obj.state_lb,[height, width, numChannel]);
                ub = reshape(obj.state_ub,[height, width, numChannel]);
                IS = ImageStar(reshape(obj.V, [height, width, numChannel,  obj.nVar + 1]), obj.C, obj.d, obj.predicate_lb, obj.predicate_ub, lb, ub);
            else
                IS = ImageStar(reshape(obj.V, [height, width, numChannel,  obj.nVar + 1]), obj.C, obj.d, obj.predicate_lb, obj.predicate_ub);
            end
                
        end

        % convert to VolumeStar set
        function VS = toVolumeStar(obj, height, width, depth, numChannel)
            % @height: height of VolumeStar
            % @width:  width of VolumeStar
            % @depth:  depth of VolumeStar
            % @numChannel: number of channels in VolumeStar
            
            if obj.dim ~= height*width*depth*numChannel
                error('Inconsistent dimension in the VolumeStar and the original Star set');
            end
            
            if ~isempty(obj.state_lb) && ~isempty(obj.state_ub)
                lb = reshape(obj.state_lb,[height, width, depth, numChannel]);
                ub = reshape(obj.state_ub,[height, width, depth, numChannel]);
                VS = VolumeStar(reshape(obj.V, [height, width, depth, numChannel,  obj.nVar + 1]), obj.C, obj.d, obj.predicate_lb, obj.predicate_ub, lb, ub);
            else
                VS = VolumeStar(reshape(obj.V, [height, width, depth, numChannel,  obj.nVar + 1]), obj.C, obj.d, obj.predicate_lb, obj.predicate_ub);
            end
                
        end

        % reset a row of a star set to zero
        function S = resetRow(obj, map)
            % @map: an array of indexes
            % xmax: max values of x[indexes]
            
            % author: Dung Tran
            % date: 6/11/2018
            
            V1 = obj.V;
            V1(map, :) = 0;
            if ~isempty(obj.Z)
                c2 = obj.Z.c;
                c2(map, :) = 0;
                V2 = obj.Z.V;
                V2(map, :) = 0;
                new_Z = Zono(c2, V2);
            else
                new_Z = [];
            end
            S = Star(V1, obj.C, obj.d, obj.predicate_lb, obj.predicate_ub, new_Z);          
        end
        
        % scale a row of a star set
        function S = scaleRow(obj, map, gamma)
            % @map: an array of indexes
            % gamma: scale value
            
            % author: Dung Tran
            % date: 11/24/2020
                
            V1 = obj.V;
            V1(map, :) = gamma*V1(map,:);
            if ~isempty(obj.Z)
                c2 = obj.Z.c;
                c2(map) = gamma*c2;
                V2 = obj.Z.V;
                V2(map, :) = gamma*V2(map, :);
                new_Z = Zono(c2, V2);
            else
                new_Z = [];
            end
            S = Star(V1, obj.C, obj.d, obj.predicate_lb, obj.predicate_ub, new_Z);          
        end
        
        % concatenate with other star
        function S = concatenate(obj, X)
            % @X: input star
            % @S: output star with is the concatenation of obj and X
            
            % author: Dung Tran
            % date: 11/18/2018
            
            if ~isa(X, 'Star')
                error('Input is not a Star');
            end
            
            c1 = obj.V(:, 1);
            c2 = X.V(:, 1);
            
            V1 = obj.V(:, 2:obj.nVar+1);
            V2 = X.V(:, 2:X.nVar + 1);
            
            c = vertcat(c1, c2);
            
            V1 = vertcat(V1, zeros(X.dim, obj.nVar));
            V2 = vertcat(zeros(obj.dim, X.nVar), V2);
            
            new_V = horzcat(c, V1, V2);
            new_C = blkdiag(obj.C, X.C);
            new_d = vertcat(obj.d, X.d);
            
            new_predicate_lb = [obj.predicate_lb; X.predicate_lb];
            new_predicate_ub = [obj.predicate_ub; X.predicate_ub];
            
            S = Star(new_V, new_C, new_d, new_predicate_lb, new_predicate_ub);            
            
        end
        
        % concatenate a star with a vector
        function S = concatenate_with_vector(obj, v)
            % @v: a vector
            % @S: output star after concatenation
            
            % author: Dung Tran
            % date: 10/1/2019
            
            
            if size(v, 2) ~= 1
                error('Input is not a vector');
            end
            
            new_c = [v; obj.V(:, 1)];
            new_V = [zeros(size(v,1), obj.nVar); obj.V(:,2:obj.nVar + 1)];          
            S = Star([new_c new_V], obj.C, obj.d, obj.predicate_lb, obj.predicate_ub);
                        
        end
        
        % change variable precision
        function S = changeVarsPrecision(obj, precision)
            S = obj;
            if strcmp(precision, 'single')
                S.V = single(S.V);
                S.C = single(S.C);
                S.d = single(S.d);
                S.predicate_lb = single(S.predicate_lb); 
                S.predicate_ub = single(S.predicate_ub);
                S.state_lb = single(S.state_lb); 
                S.state_ub = single(S.state_ub);
            elseif strcmp(precision, 'double')
                S.V = double(S.V);
                S.C = double(S.C);
                S.d = double(S.d);
                S.predicate_lb = double(S.predicate_lb); 
                S.predicate_ub = double(S.predicate_ub);
                S.state_lb = double(S.state_lb); 
                S.state_ub = double(S.state_ub);
            else
                error("Only single or double precision arrays allowed. GpuArray/dlarray are coming.")
            end
        end
        
        % change device target for the set
        function S = changeDevice(obj, deviceTarget)
            S = obj;
            if strcmp(deviceTarget, 'cpu')
                S.V = gather(S.V);
                S.C = gather(S.C);
                S.d = gather(S.d);
                S.predicate_lb = gather(S.predicate_lb); 
                S.predicate_ub = gather(S.predicate_ub);
                S.state_lb = gather(S.state_lb); 
                S.state_ub = gather(S.state_ub);
            elseif strcmp(deviceTarget, 'gpu')
                S.V = gpuArray(S.V);
                S.C = gpuArray(S.C);
                S.d = gpuArray(S.d);
                S.predicate_lb = gpuArray(S.predicate_lb); 
                S.predicate_ub = gpuArray(S.predicate_ub);
                S.state_lb = gpuArray(S.state_lb); 
                S.state_ub = gpuArray(S.state_ub);
            else
                error("Device target must be 'cpu' or 'gpu'.");
            end

        end

    end
    

    methods(Static) % plot methods

         % plot star set
        function plot(varargin)
            % author: Dung Tran
            % date: 3/27/2019
            % update: 4/2/2020
            
            switch nargin
                case 2
                    obj = varargin{1};
                    color = varargin{2};
                    map_mat = [];
                    map_vec = [];
                    approx = [];
                case 1
                    obj = varargin{1};
                    color = 'red';
                    map_mat = [];
                    map_vec = [];
                    approx = [];
                case 3
                    obj = varargin{1};
                    color = varargin{2};
                    map_mat = varargin{3};
                    map_vec = [];
                    approx = [];
                case 4
                    obj = varargin{1};
                    color = varargin{2};
                    map_mat = varargin{3};
                    map_vec = varargin{4};
                    approx = [];
                case 5
                    obj = varargin{1};
                    color = varargin{2};
                    map_mat = varargin{3};
                    map_vec = varargin{4};
                    approx = varargin{5};
                    
                otherwise
                    error('Invalid number of input arguments, should be 1, 2, 3, 4, or 5');
            end

            if isempty(map_mat)
                if obj.dim > 3
                    error('Cannot visualize the star set in > 3 dimensional space, please plot a projection of the star set');
                end
                if obj.nVar > 20
                    if isempty(approx)
                        try 
                            P = obj.toPolyhedron;
                            P.plot('color', color);
                        catch
                            warning('There was an error with MPT plotting');
                            warning('NNV plots an over-approximation of the star set using a box instead');
                            B = obj.getBox;
                            B.plot;
                        end                       
                    else                       
                        if strcmp(approx, 'zonotope')                            
                            if isempty(obj.Z)
                                Z1 = obj.getZono;
                                Zono.plots(Z1);
                            else
                                Zono.plots(obj.Z);
                            end
                        elseif strcmp(approx, 'box')
                            B = obj.getBox;
                            B.plot;
                        else
                            error('Unknown plotting option');
                        end
                    end
                else
                    if isa(obj.V, 'single') || isa(obj.C,'single') || isa(obj.d, 'single')
                        S = obj.changeVarsPrecision('double');
                        P = S.toPolyhedron;
                        P.plot('color', color);
                    else
                        P = obj.toPolyhedron;
                        P.plot('color', color);
                    end
                end
            
            else
                
                [n1,n2] = size(map_mat);
                if n1 > 3 || n1 < 1
                    error('Invalid projection matrix');
                end
                if n2 ~= obj.dim
                    error('Inconsistency between projection matrix and the star set dimension');
                end
                
                if ~isempty(map_vec)
                    [m1,m2] = size(map_vec);
                    if n1~=m1
                        error('Inconsistency between projection matrix and projection vector');
                    end
                    if m2 ~= 1
                        error('Invalid projection vector');
                    end
                end
                
                S1 = obj.affineMap(map_mat, map_vec);
                S1.plot(color,[],[],approx);
            end
            
        end
        
        % plot an array of Star (plot exactly, this is time consuming)
        function plots(varargin)
            % @S: an array of Stars
            % @color: color
            
            switch nargin
                
                case 2
                    S = varargin{1};
                    color = varargin{2};
                case 1
                    S = varargin{1};
                    color = 'b';
                otherwise
                    error('Invalid number of inputs, should be 1 or 2');
            end
            
            n = length(S);
            if n==1
                Star.plot(S,color);
            else
                for i=1:n-1
                    Star.plot(S(i),color);
                    hold on;
                end
                Star.plot(S(n), color);
            end
            
        end
        
        % plot ranges of stars at specific dimension versus time series
        function plotRanges_2D(stars, index, times, color)
            % @stars: an array of stars
            % @index: index of the specified dimension 
            % @times: an array of time points
            % @color: color of the plot
            
            % author: Dung Tran
            % date: 11/21/2018
            
            if length(stars) ~= length(times)
                error('Inconsistent length between stars and times');
            end
            
            n = length(stars);
            [ymin, ymax] = stars(1).getRange(index);
            y = [ymin ymax];
            x = [times(1) times(1)];
            plot(x, y, color);
            hold on;
            
            for i=2:n
                [ymin, ymax] = stars(i).getRange(index);
                y = [ymin ymin ymax ymax ymin];
                x = [times(i-1) times(i) times(i) times(i-1) times(i-1)];
                plot(x, y, color)
                hold on; 
            end
            hold off;
            
        end
        
        % plot an array of Stars using 2D boxes
        function plotBoxes_2D(stars, x_pos, y_pos, color)
            % plot stars using two dimension boxes
            % author: Dung Tran
            % date: 11/16/2018
                     
            n = length(stars);
            xmin = zeros(n,1);
            xmax = zeros(n,1);
            ymin = zeros(n,1);
            ymax = zeros(n,1);
            
            for i=1:n
                if isa(stars(i), 'Star')
                    [xmin(i), xmax(i)] = stars(i).getRange(x_pos);
                    [ymin(i), ymax(i)] = stars(i).getRange(y_pos);
                else
                    error('%d th input object is not a star', i);
                end
            end
                        
            for i=1:n
                x = [xmin(i) xmax(i) xmax(i) xmin(i)];
                y = [ymin(i) ymin(i) ymax(i) ymax(i)];
                hold on;
                patch(x, y, color);
                                
            end
             
        end
        
        % plot list of boxes in 2D with no fill
        function plotBoxes_2D_noFill(stars, x_pos, y_pos, color)
            % plot stars using two dimension boxes without filling
            % author: Dung Tran
            % date: 11/16/2018
            
            n = length(stars);
            xmin = zeros(n,1);
            xmax = zeros(n,1);
            ymin = zeros(n,1);
            ymax = zeros(n,1);
            
            for i=1:n
                if isa(stars(i), 'Star')
                    [xmin(i), xmax(i)] = stars(i).getRange(x_pos);
                    [ymin(i), ymax(i)] = stars(i).getRange(y_pos);
                else
                    error('%d th input object is not a star', i);
                end
            end
            
            for i=1:n
                x = [xmin(i) xmax(i) xmax(i) xmin(i) xmin(i)];
                y = [ymin(i) ymin(i) ymax(i) ymax(i) ymin(i)];
                hold on;
                plot(x, y, color);
            end
            
        end
        
        % plot list of boxes in 3D
        function plotBoxes_3D(stars, x_pos, y_pos, z_pos, color)
            % plot stars using three dimensionall boxes
            % author: Dung Tran
            % date: 11/16/2018
            
            n = length(stars);
            xmin = zeros(n,1);
            xmax = zeros(n,1);
            ymin = zeros(n,1);
            ymax = zeros(n,1);
            zmin = zeros(n,1);
            zmax = zeros(n,1);
            
            for i=1:n
                if isa(stars(i), 'Star')
                    [xmin(i), xmax(i)] = stars(i).getRange(x_pos);
                    [ymin(i), ymax(i)] = stars(i).getRange(y_pos);
                    [zmin(i), zmax(i)] = stars(i).getRange(z_pos);
                else
                    error('%d th input object is not a star', i);
                end
            end
            
            for i=1:n
                                
                p1 = [xmin(i) ymin(i) zmin(i)];
                p2 = [xmax(i) ymin(i) zmin(i)];
                p3 = [xmax(i) ymin(i) zmax(i)];
                p4 = [xmin(i) ymin(i) zmax(i)];
                p5 = [xmax(i) ymax(i) zmin(i)];
                p6 = [xmax(i) ymax(i) zmax(i)];
                p7 = [xmin(i) ymax(i) zmax(i)];
                p8 = [xmin(i) ymax(i) zmin(i)];
                
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


    methods(Static) % operate on array of stars

        % bound a set of stars by a box
        function B = get_hypercube_hull(stars)
            % @stars: an array of stars
            % @S: a box (represented as a star)
            
            n = length(stars);
            dim = stars(1).dim;
            
            for i=2:n
                if ~isa(stars(i), 'Star')
                    error('The %d th object is not a star', i);
                end                
                if stars(i).dim ~= dim
                    error('Inconsistent dimensions between stars');
                end
            end
            
            lb = [];
            ub = [];
            for i=1:n
                B1 = stars(i).getBox();
                lb = [lb B1.lb];
                ub = [ub B1.ub];
            end
            
            lb = min(lb, [], 2);
            ub = max(ub, [], 2);
            B = Box(lb, ub);            
            
        end
        
        % convex hull of stars
        function S = get_convex_hull(stars)
            % @stars: an array of stars
            % @S: a new star which is a convex hull of the input stars
            
            % author: Dung Tran
            % date: 1/21/2019
            
            n = length(stars);
            dim = stars(1).dim; 
            
            for i=2:n
                if ~isa(stars(i), 'Star')
                    error('The %d th object is not a star', i);
                end
                if stars(i).dim ~= dim
                    error('Inconsistent dimensions between stars');
                end
            end
            
            for i=1:n
                X(i) = stars(i).toPolyhedron;
            end
            
            U = PolyUnion(X);
            S = U.convexHull;
            
        end
        
        % concatenate many stars
        function S = concatenateStars(stars)
            % @stars: an array of stars
            
            new_c = [];
            new_V = [];
            new_C = [];
            new_d = [];
            new_pred_lb = [];
            new_pred_ub = [];
            
            n = length(stars);
            
            for i=1:n
                if ~isa(stars(i), 'Star')
                    error('The %d th input is not a Star', i);
                end
                
                new_c = vertcat(new_c, stars(i).V(:,1));
                new_V = blkdiag(new_V, stars(i).V(:, 2:stars(i).nVar + 1));
                new_C = blkdiag(new_C, stars(i).C);
                new_d = vertcat(new_d, stars(i).d);
                new_pred_lb = vertcat(new_pred_lb, stars(i).predicate_lb);
                new_pred_ub = vertcat(new_pred_ub, stars(i).predicate_ub);
                
            end
            
            S = Star([new_c new_V], new_C, new_d, new_pred_lb, new_pred_ub);
           
        end
       
        % merge stars using boxes and overlapness
        function S = merge_stars(I, nS, parallel)
            % @I: array of stars
            % @nP: number of stars of the output S
            
            % author: Dung Tran
            % date: 2/25/2019
            
            n = length(I);
            B = [];
            if strcmp(parallel, 'single')             
                
                for i=1:n
                    B = [B I(i).getBox];
                end

                m = I(1).dim;
                n = length(B);
                C = zeros(n, 2*m);
                for i=1:n
                    C(i, :) = [B(i).lb' B(i).ub'];
                end

                idx = kmeans(C, nS); % clustering boxes into nP groups
                R = cell(nS, 1);
                for i=1:nS
                    for j=1:n
                        if idx(j) == i
                            R{i, 1} = [R{i, 1} B(j)];
                        end
                    end
                end

                S = [];
                for i=1:nS
                    B = Box.boxHull(R{i, 1}); % return a box                    
                    S = [S B.toStar];
                end
                
            elseif strcmp(parallel, 'parallel')
                
                parfor i=1:n
                    B = [B I(i).getBox];
                end
                m = I(1).dim;
                n = length(B);
                C = zeros(n, 2*m);
                for i=1:n
                    C(i, :) = [B(i).lb' B(i).ub'];
                end
                idx = kmeans(C, nS); % clustering boxes into nP groups
                R = cell(nS, 1);
                for i=1:nS
                    for j=1:n
                        if idx(j) == i
                            R{i, 1} = [R{i, 1} B(j)];
                        end
                    end
                end

                S = [];
                parfor i=1:nS
                    B = Box.boxHull(R{i, 1}); % return a box                    
                    S = [S B.toStar];
                end
                
            else
                error('Unknown parallel option');
            end

        end
        
    end


    methods(Static) % helper functions

        % generate random star set
        function S = rand(dim)
            % @dim: dimension of the random star set
            % @S: the star set
            
            % author: Dung Tran
            % date: 9/16/2020
            
            if dim <= 0 
                error('Invalid dimension');
            end
            P = ExamplePoly.randHrep('d',dim); % random polyhedron
            S = Star(P);  
            
        end
           
    end
    
end