classdef GCNLayer < handle
    % The GCNLayer class for Graph Convolutional Networks
    %   Performs: Y = A_norm * X * W + b
    %   Where:
    %     A_norm: Normalized adjacency matrix (N x N)
    %     X: Input node features (N x F_in)
    %     W: Weight matrix (F_in x F_out)
    %     b: Bias vector (F_out x 1)
    %     Y: Output node features (N x F_out)
    %
    %   Note: ReLU activation is applied by a separate ReluLayer
    %
    %   Reference: Kipf & Welling, "Semi-Supervised Classification with
    %              Graph Convolutional Networks", ICLR 2017
    %
    %   Anne Tumlin: 1/13/2026

    properties
        Name = 'gcn_layer';
        InputSize = 0;   % F_in: number of input features per node
        OutputSize = 0;  % F_out: number of output features per node
        Weights = [];    % Weight matrix (F_in x F_out)
        Bias = [];       % Bias vector (F_out x 1)

        % Standard layer interface
        NumInputs = 1;
        InputNames = {'in'};
        NumOutputs = 1;
        OutputNames = {'out'};
    end


    methods % main methods

        % constructor of the class
        function obj = GCNLayer(varargin)
            % GCNLayer constructor
            % Usage:
            %   GCNLayer() - empty layer
            %   GCNLayer(W, b) - with weights and bias
            %   GCNLayer(name, W, b) - with name, weights, and bias
            %
            % @W: Weight matrix (F_in x F_out)
            % @b: Bias vector (F_out x 1)

            % author: Anne Tumlin
            % date: 1/13/2026

            switch nargin

                case 3
                    name = varargin{1};
                    W = varargin{2};
                    b = varargin{3};
                    if ~ischar(name) && ~isstring(name)
                        error('Name must be a string or char array');
                    else
                        obj.Name = name;
                    end
                    if size(W, 2) ~= size(b, 1)
                        error('Inconsistent dimensions: size(W,2) must equal size(b,1)');
                    end
                    if size(b, 2) ~= 1
                        error('Bias vector should have one column');
                    end

                    obj.InputSize = size(W, 1);
                    obj.OutputSize = size(W, 2);
                    obj.Weights = W;
                    obj.Bias = b;

                case 2

                    W = varargin{1};
                    b = varargin{2};
                    obj.Name = 'gcn_layer';
                    if size(W, 2) ~= size(b, 1)
                        error('Inconsistent dimensions: size(W,2) must equal size(b,1)');
                    end
                    if size(b, 2) ~= 1
                        error('Bias vector should have one column');
                    end
                    obj.InputSize = size(W, 1);
                    obj.OutputSize = size(W, 2);
                    obj.Weights = W;
                    obj.Bias = b;

                case 0

                    obj.Name = 'gcn_layer';
                    obj.InputSize = 0;
                    obj.OutputSize = 0;
                    obj.Weights = [];
                    obj.Bias = [];

                otherwise
                    error('Invalid number of inputs (should be 0, 2, or 3)');
            end

        end

        % evaluation method
        function Y = evaluate(obj, X, A_norm)
            % Forward pass through GCN layer
            % @X: Input node features (N x F_in)
            % @A_norm: Normalized adjacency matrix (N x N)
            % @Y: Output node features (N x F_out)

            % author: Anne Tumlin
            % date: 1/13/2026

            if size(X, 2) ~= obj.InputSize
                error('Input feature dimension mismatch: expected %d, got %d', obj.InputSize, size(X, 2));
            end

            if size(A_norm, 1) ~= size(X, 1) || size(A_norm, 2) ~= size(X, 1)
                error('Adjacency matrix dimension mismatch with number of nodes');
            end

            % Graph convolution: aggregate then transform
            X_agg = A_norm * X;              % Aggregate neighbors: N x F_in
            Y = X_agg * obj.Weights;         % Linear transform: N x F_out

            if ~isempty(obj.Bias)
                Y = Y + obj.Bias';           % Broadcast bias across nodes
            end

        end

        % main reachability analysis function
        function S = reach(varargin)
            % Reachability analysis for GCNLayer
            % @in_set: Input GraphStar set
            % @A_norm: Normalized adjacency matrix
            % @method: 'approx-star', 'exact-star', 'approx-zono'
            % @option: 'single' or 'parallel'

            % author: Anne Tumlin
            % date: 1/13/2026

            switch nargin

                case 7
                    obj = varargin{1};
                    in_sets = varargin{2};
                    A_norm = varargin{3};
                    method = varargin{4};
                    option = varargin{5};
                    % relaxFactor = varargin{6}; % not used
                    % lp_solver = varargin{7}; % not used

                case 6
                    obj = varargin{1};
                    in_sets = varargin{2};
                    A_norm = varargin{3};
                    method = varargin{4};
                    option = varargin{5};
                    % relaxFactor = varargin{6}; % not used

                case 5
                    obj = varargin{1};
                    in_sets = varargin{2};
                    A_norm = varargin{3};
                    method = varargin{4};
                    option = varargin{5};

                case 4
                    obj = varargin{1};
                    in_sets = varargin{2};
                    A_norm = varargin{3};
                    method = varargin{4};
                    option = [];

                case 3
                    obj = varargin{1};
                    in_sets = varargin{2};
                    A_norm = varargin{3};
                    method = 'approx-star';
                    option = [];

                otherwise
                    error('Invalid number of input arguments (should be 2-6)');
            end

            if strcmp(method, 'approx-star') || strcmp(method, 'exact-star') || strcmp(method, 'abs-dom') || contains(method, "relax-star")
                S = obj.reach_star_multipleInputs(in_sets, A_norm, option);
            elseif strcmp(method, 'approx-zono')
                S = obj.reach_zono_multipleInputs(in_sets, A_norm, option);
            else
                error('Unknown reachability method: %s', method);
            end

        end

    end


    methods % reachability methods

        % reachability for a single GraphStar input
        function gs_out = reach_star_single_input(obj, in_gs, A_norm)
            % Reachability through GCN layer for single GraphStar
            % @in_gs: Input GraphStar set
            % @A_norm: Normalized adjacency matrix (N x N)
            % @gs_out: Output GraphStar set

            % author: Anne Tumlin
            % date: 1/13/2026

            if ~isa(in_gs, 'GraphStar')
                error('Input must be a GraphStar');
            end

            if in_gs.numFeatures ~= obj.InputSize
                error('Input feature dimension mismatch: expected %d, got %d', obj.InputSize, in_gs.numFeatures);
            end

            numNodes = in_gs.numNodes;
            numPred = in_gs.numPred;

            % Preallocate output V matrix
            V_out = zeros(numNodes, obj.OutputSize, numPred + 1, 'like', in_gs.V);

            % Transform each generator through GCN operation
            for k = 1:(numPred + 1)
                X_k = in_gs.V(:, :, k);           % [N x F_in]
                X_agg = A_norm * X_k;              % Aggregate neighbors
                X_trans = X_agg * obj.Weights;    % Linear transform

                if k == 1 && ~isempty(obj.Bias)
                    X_trans = X_trans + obj.Bias'; % Bias only for center
                end

                V_out(:, :, k) = X_trans;
            end

            % Constraints unchanged (linear transformation)
            gs_out = GraphStar(V_out, in_gs.C, in_gs.d, in_gs.pred_lb, in_gs.pred_ub);

        end

        % reachability for multiple GraphStar inputs
        function S = reach_star_multipleInputs(obj, in_sets, A_norm, option)
            % Reachability for multiple GraphStar inputs
            % @in_sets: Array of GraphStar sets
            % @A_norm: Normalized adjacency matrix
            % @option: 'single' or 'parallel'

            % author: Anne Tumlin
            % date: 1/13/2026

            n = length(in_sets);

            if n == 1
                S = obj.reach_star_single_input(in_sets, A_norm);
            else
                S(n) = GraphStar;

                if strcmp(option, 'parallel')
                    parfor i = 1:n
                        S(i) = obj.reach_star_single_input(in_sets(i), A_norm);
                    end
                else
                    for i = 1:n
                        S(i) = obj.reach_star_single_input(in_sets(i), A_norm);
                    end
                end
            end

        end

        % % zonotope reachability (placeholder)
        % function S = reach_zono_multipleInputs(obj, in_sets, A_norm, option)
        %     % Zonotope-based reachability (not yet implemented)

        %     error('Zonotope reachability for GCNLayer not yet implemented');

        % end

    end


    methods % helper methods

        function obj = toGPU(obj)
            % Move weights to GPU

            obj.Weights = gpuArray(obj.Weights);
            if ~isempty(obj.Bias)
                obj.Bias = gpuArray(obj.Bias);
            end
        end

        function obj = toCPU(obj)
            % Move weights to CPU

            obj.Weights = gather(obj.Weights);
            if ~isempty(obj.Bias)
                obj.Bias = gather(obj.Bias);
            end
        end

        function obj = changeParamsPrecision(obj, precision)
            % Change precision of parameters
            % @precision: 'single' or 'double'

            if strcmp(precision, 'double')
                obj.Weights = double(obj.Weights);
                if ~isempty(obj.Bias)
                    obj.Bias = double(obj.Bias);
                end
            elseif strcmp(precision, 'single')
                obj.Weights = single(obj.Weights);
                if ~isempty(obj.Bias)
                    obj.Bias = single(obj.Bias);
                end
            else
                error('Precision must be "single" or "double"');
            end
        end

    end

end
