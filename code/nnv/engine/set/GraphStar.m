classdef GraphStar < handle
    % Class for representing set of graph node features using Star set
    % Graph node features can be perturbed by bounded noise. A perturbed
    % graph can be represented using a GraphStar Set.
    % Anne Tumlin: 1/5/2026

    %=================================================================%
    %   A graph with N nodes and F features per node is represented by
    %   an N x F matrix of node features.
    %
    %   Problem: How to represent a perturbed graph?
    %
    %   Use center node features (a matrix) + a perturbation matrix
    %   (bounds on perturbations to each node feature)
    %
    %   For example: Consider a graph with 3 nodes and 2 features
    %   The node features are represented by 3 x 2 matrix:
    %               NF = [1.0 0.5; 0.8 0.3; 0.6 0.9]
    %   This graph is perturbed at all node features by bounded noise:
    %               LB = [-0.1 -0.1; -0.1 -0.1; -0.1 -0.1]
    %               UB = [0.1 0.1; 0.1 0.1; 0.1 0.1]
    %
    %   Under perturbation: NF(i,j) + LB(i,j) <= x(i,j) <= NF(i,j) + UB(i,j)
    %
    %   To represent the perturbed graph we use NF, LB, UB matrices
    %   The graph object is: graph = GraphStar(NF, LB, UB)
    %=================================================================%


    properties
        numNodes = 0;      % number of nodes in graph
        numFeatures = 0;   % number of features per node

        % A box representation of a GraphStar
        % A convenient way for user to specify the perturbation

        NF = [];  % center node features (N x F matrix)
        LB = [];  % lower bound of perturbation (N x F matrix)
        UB = [];  % upper bound of perturbation (N x F matrix)

        % Star representation of a GraphStar
        % ====================================================================%
        %                   Definition of GraphStar
        %
        % A GraphStar set S is defined by:
        % S = {x| x = V[0] + a[1]*V[1] + a[2]*V[2] + ... + a[n]*V[n]
        %           = V * b, V = {c V[1] V[2] ... V[n]},
        %                    b = [1 a[1] a[2] ... a[n]]^T
        %                    where C*a <= d, constraints on a[i]}
        % where, V[0], V[i] are 2D matrices with the same dimension, i.e.,
        % V[i] \in R^{N x F}
        % V[0] : is called the center matrix and V[i] is called the basis matrix
        % [a[1]...a[n] are called predicate variables
        % C: is the predicate constraint matrix
        % d: is the predicate constraint vector
        %
        % ====================================================================%
        % The Star representation is convenient for reachability analysis
        V = [];        % basis matrices [numNodes x numFeatures x (numPred+1)]
        C = [];        % constraint matrix of the predicates
        d = [];        % constraint vector of the predicates
        numPred = 0;   % number of predicate variables
        pred_lb = [];  % lower bound vector of the predicates
        pred_ub = [];  % upper bound vector of the predicates
        nf_lb = [];    % lower bound node features of the GraphStar
        nf_ub = [];    % upper bound node features of the GraphStar

    end

    methods % Constructor and sampling methods

        % constructor using box representation or Star representation
        function obj = GraphStar(varargin)
            % @nargin = 3: NF = varargin{1}, LB = varagin{2}, UB = varargin{3}
            %         = 2: lb = varargin{1}, ub = varargin{2}
            %         = 5: V, C, d, pred_lb, pred_ub
            %         = 7: V, C, d, pred_lb, pred_ub, nf_lb, nf_ub
            %         = 0: empty GraphStar

            % author: Anne Tumlin
            % date: 1/5/2026

            switch nargin

                case 3 % input center node features and perturbation bound matrices

                    NF1 = varargin{1};
                    LB1 = varargin{2};
                    UB1 = varargin{3};
                    n = size(NF1);
                    l = size(LB1);
                    u = size(UB1);

                    if length(n) ~= 2 || length(l) ~= 2 || length(u) ~= 2
                        error('Node feature matrices must be 2D (numNodes x numFeatures)');
                    end

                    if n(1) ~= l(1) || n(1) ~= u(1) || n(2) ~= l(2) || n(2) ~= u(2)
                        error('Inconsistency between center node features and perturbation bound matrices');
                    end

                    obj.numNodes = n(1);
                    obj.numFeatures = n(2);
                    obj.NF = NF1;
                    obj.LB = LB1;
                    obj.UB = UB1;

                    obj.nf_lb = NF1 + LB1;  % lower bound node features
                    obj.nf_ub = NF1 + UB1;  % upper bound node features

                    % converting box GraphStar to Star representation
                    N = obj.numNodes * obj.numFeatures;
                    I = Star(reshape(obj.nf_lb, [N, 1]), reshape(obj.nf_ub, [N, 1]));
                    obj.V = reshape(I.V, [obj.numNodes, obj.numFeatures, I.nVar + 1]);
                    obj.C = I.C;
                    obj.d = I.d;
                    obj.pred_lb = I.predicate_lb;
                    obj.pred_ub = I.predicate_ub;
                    obj.numPred = I.nVar;

                case 5 % Star representation

                    V1 = varargin{1};      % basis matrices
                    C1 = varargin{2};      % predicate constraint matrix
                    d1 = varargin{3};      % predicate constraint vector
                    pred_lb1 = varargin{4}; % predicate lower bound
                    pred_ub1 = varargin{5}; % predicate upper bound

                    if size(C1, 1) ~= size(d1, 1)
                        error('Inconsistent dimension between constraint matrix and constraint vector');
                    end

                    if size(d1, 2) ~= 1
                        error('Invalid constraint vector, vector should have one column');
                    end

                    obj.numPred = size(C1, 2);
                    obj.C = C1;
                    obj.d = d1;

                    if size(C1, 2) ~= size(pred_lb1, 1) || size(C1, 2) ~= size(pred_ub1, 1)
                        error('Number of predicates is different from the size of the lower bound or upper bound predicate vector');
                    end

                    if size(pred_lb1, 2) ~= 1 || size(pred_ub1, 2) ~= 1
                        error('Invalid lower/upper bound predicate vector, vector should have one column');
                    end

                    obj.pred_lb = pred_lb1;
                    obj.pred_ub = pred_ub1;

                    n = size(V1);

                    if length(n) == 2
                        % V1 is [numNodes x numFeatures], single center (no predicates case)
                        obj.V = V1;
                        obj.numNodes = n(1);
                        obj.numFeatures = n(2);
                    elseif length(n) == 3
                        if n(3) ~= obj.numPred + 1
                            error('Inconsistency between the basis matrix and the number of predicate variables');
                        else
                            obj.V = V1;
                            obj.numNodes = n(1);
                            obj.numFeatures = n(2);
                        end
                    else
                        error('Invalid basis matrix, should be 2D or 3D');
                    end

                case 7 % Star representation with cached bounds

                    V1 = varargin{1};       % basis matrices
                    C1 = varargin{2};       % predicate constraint matrix
                    d1 = varargin{3};       % predicate constraint vector
                    pred_lb1 = varargin{4}; % predicate lower bound
                    pred_ub1 = varargin{5}; % predicate upper bound
                    nf_lb1 = varargin{6};   % lower bound node features
                    nf_ub1 = varargin{7};   % upper bound node features

                    if size(C1, 1) ~= size(d1, 1)
                        error('Inconsistent dimension between constraint matrix and constraint vector');
                    end

                    if size(d1, 2) ~= 1
                        error('Invalid constraint vector, vector should have one column');
                    end

                    obj.numPred = size(C1, 2);
                    obj.C = C1;
                    obj.d = d1;

                    if size(C1, 2) ~= size(pred_lb1, 1) || size(C1, 2) ~= size(pred_ub1, 1)
                        error('Number of predicates is different from the size of the lower bound or upper bound predicate vector');
                    end

                    if size(pred_lb1, 2) ~= 1 || size(pred_ub1, 2) ~= 1
                        error('Invalid lower/upper bound predicate vector, vector should have one column');
                    end

                    obj.pred_lb = pred_lb1;
                    obj.pred_ub = pred_ub1;

                    n = size(V1);

                    if length(n) == 2
                        obj.V = V1;
                        obj.numNodes = n(1);
                        obj.numFeatures = n(2);
                    elseif length(n) == 3
                        if n(3) ~= obj.numPred + 1
                            error('Inconsistency between the basis matrix and the number of predicate variables');
                        else
                            obj.V = V1;
                            obj.numNodes = n(1);
                            obj.numFeatures = n(2);
                        end
                    else
                        error('Invalid basis matrix, should be 2D or 3D');
                    end

                    obj.nf_lb = nf_lb1;
                    obj.nf_ub = nf_ub1;

                case 2 % create GraphStar from lower and upper bounds directly

                    lb_nf = varargin{1};
                    ub_nf = varargin{2};

                    n = size(lb_nf);
                    m = size(ub_nf);

                    if length(n) ~= 2 || length(m) ~= 2
                        error('Bound matrices must be 2D (numNodes x numFeatures)');
                    end

                    if n(1) ~= m(1) || n(2) ~= m(2)
                        error('Inconsistency between lower bound and upper bound matrices');
                    end

                    obj.numNodes = n(1);
                    obj.numFeatures = n(2);

                    N = n(1) * n(2);
                    lb = reshape(lb_nf, [N, 1]);
                    ub = reshape(ub_nf, [N, 1]);
                    S = Star(lb, ub);

                    obj.V = reshape(S.V, [n(1), n(2), S.nVar + 1]);
                    obj.C = S.C;
                    obj.d = S.d;
                    obj.pred_lb = S.predicate_lb;
                    obj.pred_ub = S.predicate_ub;
                    obj.numPred = S.nVar;
                    obj.nf_lb = lb_nf;
                    obj.nf_ub = ub_nf;

                case 0 % create an empty GraphStar

                    obj.numNodes = 0;
                    obj.numFeatures = 0;
                    obj.NF = [];
                    obj.LB = [];
                    obj.UB = [];
                    obj.V = [];
                    obj.C = [];
                    obj.d = [];
                    obj.numPred = 0;
                    obj.pred_lb = [];
                    obj.pred_ub = [];
                    obj.nf_lb = [];
                    obj.nf_ub = [];

                otherwise
                    error('Invalid number of input arguments, (should be 0, 2, 3, 5, or 7)');
            end

        end

        % randomly generate a set of node feature matrices from a GraphStar
        function node_features = sample(obj, N)
            % @N: number of samples

            % author: Anne Tumlin
            % date: 1/5/2026

            if isempty(obj.V)
                error('The GraphStar is an empty set');
            end

            if isempty(obj.C) || isempty(obj.d)
                node_features = {obj.NF};
            else
                V1 = [cast(zeros(obj.numPred, 1), 'like', obj.C)  eye(obj.numPred)];
                S = Star(V1, obj.C, obj.d, obj.pred_lb, obj.pred_ub);
                pred_samples = S.sample(N);

                M = size(pred_samples, 2);
                node_features = cell(1, M);
                for i=1:M
                    node_features{i} = obj.evaluate(pred_samples(:, i));
                end
            end

        end

        % evaluate a GraphStar with specific values of predicates
        function nf = evaluate(obj, pred_val)
            % @pred_val: valued vector of predicate variables

            % author: Anne Tumlin
            % date: 1/5/2026

            if isempty(obj.V)
                error('The GraphStar is an empty set');
            end

            if size(pred_val, 2) ~= 1
                error('Invalid predicate vector');
            end

            if size(pred_val, 1) ~= obj.numPred
                error('Inconsistency between the size of the predicate vector and the number of predicates in the GraphStar');
            end

            nf = obj.V(:, :, 1);
            for j=2:obj.numPred + 1
                nf = nf + pred_val(j-1) * obj.V(:, :, j);
            end

        end

    end


    methods % operations / transformation functions

        % affineMap of a GraphStar is another GraphStar (y = scale * x + offset)
        function gs = affineMap(obj, scale, offset)
            % @scale: scale coefficient (scalar or matrix)
            % @offset: offset coefficient (scalar or matrix)
            % @gs: a new GraphStar

            % author: Anne Tumlin
            % date: 1/5/2026

            if ~isempty(scale)
                new_V = scale .* obj.V;
            else
                new_V = obj.V;
            end

            if ~isempty(offset)
                new_V(:, :, 1) = new_V(:, :, 1) + offset;
            end
            gs = GraphStar(new_V, obj.C, obj.d, obj.pred_lb, obj.pred_ub);

        end

        % Minkowski Sum of two GraphStars is another GraphStar (y = x1 + x2)
        function gs = MinkowskiSum(obj, GS)

            % author: Anne Tumlin
            % date: 1/5/2026

            S1 = obj.toStar;
            S2 = GS.toStar;
            S = S1.MinkowskiSum(S2);

            % Reshape back to GraphStar
            new_V = reshape(S.V, [obj.numNodes, obj.numFeatures, S.nVar + 1]);
            gs = GraphStar(new_V, S.C, S.d, S.predicate_lb, S.predicate_ub);

        end

        % transform to Star
        function S = toStar(obj)

            % author: Anne Tumlin
            % date: 1/5/2026

            N = obj.numNodes * obj.numFeatures;
            np = obj.numPred;

            V1 = cast(zeros(N, np+1), 'like', obj.V);
            for j=1:np+1
                V1(:, j) = reshape(obj.V(:, :, j), N, 1);
            end
            if ~isempty(obj.nf_lb) && ~isempty(obj.nf_ub)
                state_lb = reshape(obj.nf_lb, N, 1);
                state_ub = reshape(obj.nf_ub, N, 1);
                S = Star(V1, obj.C, obj.d, obj.pred_lb, obj.pred_ub, state_lb, state_ub);
            else
                S = Star(V1, obj.C, obj.d, obj.pred_lb, obj.pred_ub);
            end

        end

        % checking if a GraphStar is an empty set
        function bool = isEmptySet(obj)

            % author: Anne Tumlin
            % date: 1/5/2026

            try
                S = obj.toStar;
                bool = S.isEmptySet;
            catch
                [lb, ub] = obj.estimateRanges;
                if isempty(lb) && isempty(ub)
                    bool = 1;
                else
                    bool = 0;
                end
            end
        end

        % check if a GraphStar contains a specific node feature matrix
        function bool = contains(obj, nf)
            % @nf: input node features matrix
            % @bool: = 1 if the GraphStar contains nf
            %        = 0 if the GraphStar does not contain nf

            % author: Anne Tumlin
            % date: 1/5/2026

            n = size(nf);
            if length(n) ~= 2
                error('Input must be a 2D matrix');
            end

            if n(1) ~= obj.numNodes || n(2) ~= obj.numFeatures
                error('Inconsistent dimensions between input and the GraphStar');
            end

            y = reshape(nf, [n(1)*n(2), 1]);
            S = obj.toStar;
            bool = S.contains(y);

        end

    end


    methods % get methods

        % get ranges of a node feature at specific position
        function [xmin, xmax] = getRange(varargin)
            % @node_ind: node index
            % @feat_ind: feature index
            % @xmin: min of x(node_ind, feat_ind)
            % @xmax: max of x(node_ind, feat_ind)

            % author: Anne Tumlin
            % date: 1/5/2026

            switch nargin
                case 3
                    obj = varargin{1};
                    node_ind = varargin{2};
                    feat_ind = varargin{3};
                    lp_solver = 'linprog';
                case 4
                    obj = varargin{1};
                    node_ind = varargin{2};
                    feat_ind = varargin{3};
                    lp_solver = varargin{4};
                otherwise
                    error('Invalid number of input arguments, should be 2 or 3');
            end

            if isempty(obj.C) || isempty(obj.d)
                error('The GraphStar is empty');
            end

            if node_ind < 1 || node_ind > obj.numNodes
                error('Invalid node index');
            end

            if feat_ind < 1 || feat_ind > obj.numFeatures
                error('Invalid feature index');
            end

            % min
            f = reshape(obj.V(node_ind, feat_ind, 2:obj.numPred + 1), [1, obj.numPred]);
            [fval, exitflag] = lpsolver(f, obj.C, obj.d, [], [], obj.pred_lb, obj.pred_ub, lp_solver);
            if ismember(exitflag, ["l1","g5"])
               xmin = fval + obj.V(node_ind, feat_ind, 1);
            else
                error("Cannot find an optimal solution, exitflag = " + string(exitflag));
            end
            % max
            [fval, exitflag] = lpsolver(-f, obj.C, obj.d, [], [], obj.pred_lb, obj.pred_ub, lp_solver);
            if ismember(exitflag, ["l1","g5"])
                xmax = -fval + obj.V(node_ind, feat_ind, 1);
            else
                error("Cannot find an optimal solution, exitflag = " + string(exitflag));
            end

            obj.nf_lb(node_ind, feat_ind) = xmin;
            obj.nf_ub(node_ind, feat_ind) = xmax;

        end

        % estimate ranges quickly using only predicate bound information
        function [nf_lb_out, nf_ub_out] = estimateRanges(varargin)
            % @nf_lb_out: lower bound node features
            % @nf_ub_out: upper bound node features

            % author: Anne Tumlin
            % date: 1/5/2026

            switch nargin
                case 1
                    obj = varargin{1};
                otherwise
                    error('Invalid number of input arguments, should be 0');
            end

            if isempty(obj.C) || isempty(obj.d)
                warning('The GraphStar is empty');
                obj.nf_lb = obj.V(:,:,1);
                obj.nf_ub = obj.V(:,:,1);
            end

            if isempty(obj.nf_lb) || isempty(obj.nf_ub)

                gens = obj.V(:, :, 2:end);
                pos_gens = gens;
                pos_gens(gens < 0) = 0;
                neg_gens = gens;
                neg_gens(gens > 0) = 0;
                nf_lb_out = obj.V(:, :, 1) + tensorprod(pos_gens, obj.pred_lb, 3, 1) + tensorprod(neg_gens, obj.pred_ub, 3, 1);
                nf_ub_out = obj.V(:, :, 1) + tensorprod(pos_gens, obj.pred_ub, 3, 1) + tensorprod(neg_gens, obj.pred_lb, 3, 1);

                obj.nf_lb = nf_lb_out;
                obj.nf_ub = nf_ub_out;

            else

                nf_lb_out = obj.nf_lb;
                nf_ub_out = obj.nf_ub;

            end

        end

        % get lower bound and upper bound node features using LP
        function [nf_lb_out, nf_ub_out] = getRanges(varargin)
            % @nf_lb_out: lower bound node features
            % @nf_ub_out: upper bound node features

            % author: Anne Tumlin
            % date: 1/5/2026

            switch nargin
                case 1
                    obj = varargin{1};
                    lp_solver = 'linprog';
                case 2
                    obj = varargin{1};
                    lp_solver = varargin{2};
                otherwise
                    error('Invalid number of input arguments');
            end

            nf_lb_out = zeros(obj.numNodes, obj.numFeatures);
            nf_ub_out = zeros(obj.numNodes, obj.numFeatures);

            for i=1:obj.numNodes
                for j=1:obj.numFeatures
                    [nf_lb_out(i, j), nf_ub_out(i, j)] = obj.getRange(i, j, lp_solver);
                end
            end

            obj.nf_lb = nf_lb_out;
            obj.nf_ub = nf_ub_out;

        end

        % quickly estimate range for a specific node feature
        function [xmin, xmax] = estimateRange(obj, node_ind, feat_ind)
            % @node_ind: node index
            % @feat_ind: feature index
            % @xmin: min of x(node_ind, feat_ind)
            % @xmax: max of x(node_ind, feat_ind)

            % author: Anne Tumlin
            % date: 1/5/2026

            if isempty(obj.C) || isempty(obj.d)
                error('The GraphStar is empty');
            end

            if node_ind < 1 || node_ind > obj.numNodes
                error('Invalid node index');
            end

            if feat_ind < 1 || feat_ind > obj.numFeatures
                error('Invalid feature index');
            end

            f = obj.V(node_ind, feat_ind, 1:obj.numPred + 1);
            xmin = f(1);
            xmax = f(1);

            for i=2:obj.numPred + 1
                if f(i) >= 0
                    xmin = xmin + f(i) * obj.pred_lb(i-1);
                    xmax = xmax + f(i) * obj.pred_ub(i-1);
                else
                    xmin = xmin + f(i) * obj.pred_ub(i-1);
                    xmax = xmax + f(i) * obj.pred_lb(i-1);
                end
            end

        end

    end

end
