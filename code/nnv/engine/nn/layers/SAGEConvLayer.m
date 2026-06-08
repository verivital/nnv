classdef SAGEConvLayer < handle
    % The SAGEConvLayer class for GraphSAGE networks
    %   Performs: Y = X * W_node + A * X * W_edge + b
    %   Where:
    %     A: Binary adjacency matrix (N x N), NOT normalized
    %     X: Input node features (N x F_in)
    %     W_node: Node (self) weight matrix (F_in x F_out)
    %     W_edge: Edge (neighbor) weight matrix (F_in x F_out)
    %     b: Bias vector (F_out x 1)
    %     Y: Output node features (N x F_out)
    %
    %   Note: ReLU activation is applied by a separate ReluLayer
    %
    %   Weight Format: W_node and W_edge are (F_in x F_out). If importing
    %   from PyTorch or SCIP-MPNN, transpose since they use (F_out x F_in).
    %
    %   Adjacency: A is a raw binary adjacency matrix (no self-loops, no
    %   degree normalization). This matches SCIP-MPNN's SAGEConv with
    %   SUM aggregation.
    %
    % Main references:
    % 1) Hamilton et al., "Inductive Representation Learning on Large
    %    Graphs", NeurIPS 2017
    % 2) Hojny et al., "Verifying message-passing neural networks via
    %    topology-based bounds tightening", ICML 2024
    %
    % Author: Anne Tumlin
    % Date: 02/12/2026

    properties
        Name = 'sage_conv_layer';
        InputSize = 0;       % F_in: number of input features per node
        OutputSize = 0;      % F_out: number of output features per node
        NodeWeights = [];    % Node (self) weight matrix (F_in x F_out)
        EdgeWeights = [];    % Edge (neighbor) weight matrix (F_in x F_out)
        Bias = [];           % Bias vector (F_out x 1)

        % Standard layer interface
        NumInputs = 1;
        InputNames = {'in'};
        NumOutputs = 1;
        OutputNames = {'out'};
    end


    methods % main methods

        function obj = SAGEConvLayer(varargin)
            % SAGEConvLayer constructor
            % Usage:
            %   SAGEConvLayer() - empty layer
            %   SAGEConvLayer(W_node, W_edge, b) - with weights and bias
            %   SAGEConvLayer(name, W_node, W_edge, b) - with name
            %
            % @W_node: Node weight matrix (F_in x F_out)
            % @W_edge: Edge weight matrix (F_in x F_out)
            % @b: Bias vector (F_out x 1)

            switch nargin

                case 4
                    name = varargin{1};
                    W_node = varargin{2};
                    W_edge = varargin{3};
                    b = varargin{4};
                    if ~ischar(name) && ~isstring(name)
                        error('Name must be a string or char array');
                    end
                    obj.Name = name;
                    obj = obj.setWeights(W_node, W_edge, b);

                case 3
                    W_node = varargin{1};
                    W_edge = varargin{2};
                    b = varargin{3};
                    obj = obj.setWeights(W_node, W_edge, b);

                case 0
                    % Empty constructor

                otherwise
                    error('Invalid number of inputs (should be 0, 3, or 4)');
            end
        end


        function Y = evaluate(obj, X, A)
            % Forward pass through SAGEConv layer
            % @X: Input node features (N x F_in)
            % @A: Binary adjacency matrix (N x N)
            % @Y: Output node features (N x F_out)

            if size(X, 2) ~= obj.InputSize
                error('Input feature dimension mismatch: expected %d, got %d', obj.InputSize, size(X, 2));
            end

            if size(A, 1) ~= size(X, 1) || size(A, 2) ~= size(X, 1)
                error('Adjacency matrix dimension mismatch with number of nodes');
            end

            % Self-node transformation
            Y_node = X * obj.NodeWeights;       % N x F_out

            % Neighbor aggregation (SUM) then transformation
            X_agg = A * X;                      % N x F_in (sum of neighbors)
            Y_edge = X_agg * obj.EdgeWeights;   % N x F_out

            % Combine
            Y = Y_node + Y_edge;

            if ~isempty(obj.Bias)
                Y = Y + obj.Bias';              % Broadcast bias across nodes
            end
        end


        function S = reach(varargin)
            % Reachability analysis for SAGEConvLayer
            % @in_set: Input GraphStar set
            % @A: Binary adjacency matrix
            % @method: 'approx-star', 'exact-star', 'approx-zono'
            % @option: 'single' or 'parallel'

            switch nargin
                case 7
                    obj = varargin{1};
                    in_sets = varargin{2};
                    A = varargin{3};
                    method = varargin{4};
                    option = varargin{5};

                case 6
                    obj = varargin{1};
                    in_sets = varargin{2};
                    A = varargin{3};
                    method = varargin{4};
                    option = varargin{5};

                case 5
                    obj = varargin{1};
                    in_sets = varargin{2};
                    A = varargin{3};
                    method = varargin{4};
                    option = varargin{5};

                case 4
                    obj = varargin{1};
                    in_sets = varargin{2};
                    A = varargin{3};
                    method = varargin{4};
                    option = [];

                case 3
                    obj = varargin{1};
                    in_sets = varargin{2};
                    A = varargin{3};
                    method = 'approx-star';
                    option = [];

                otherwise
                    error('Invalid number of input arguments (should be 2-6)');
            end

            if strcmp(method, 'approx-star') || strcmp(method, 'exact-star') || strcmp(method, 'abs-dom') || contains(method, "relax-star")
                S = obj.reach_star_multipleInputs(in_sets, A, option);
            elseif strcmp(method, 'approx-zono')
                S = obj.reach_zono_multipleInputs(in_sets, A, option);
            else
                error('Unknown reachability method: %s', method);
            end
        end

    end


    methods % reachability methods

        function gs_out = reach_star_single_input(obj, in_gs, A)
            % Reachability through SAGEConv layer for single GraphStar
            % @in_gs: Input GraphStar set
            % @A: Binary adjacency matrix (N x N)
            % @gs_out: Output GraphStar set
            %
            % SAGEConv is linear in X for fixed A, so reachability is exact.

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

            % Transform each generator through SAGEConv operation
            for k = 1:(numPred + 1)
                X_k = in_gs.V(:, :, k);                    % [N x F_in]
                Y_node = X_k * obj.NodeWeights;             % Self-node: N x F_out
                Y_edge = (A * X_k) * obj.EdgeWeights;       % Neighbor: N x F_out
                X_trans = Y_node + Y_edge;

                if k == 1 && ~isempty(obj.Bias)
                    X_trans = X_trans + obj.Bias';           % Bias only for center
                end

                V_out(:, :, k) = X_trans;
            end

            % Constraints unchanged (linear transformation)
            gs_out = GraphStar(V_out, in_gs.C, in_gs.d, in_gs.pred_lb, in_gs.pred_ub);
        end


        function S = reach_star_multipleInputs(obj, in_sets, A, option)
            % Reachability for multiple GraphStar inputs

            n = length(in_sets);

            if n == 1
                S = obj.reach_star_single_input(in_sets, A);
            else
                S(n) = GraphStar;

                if strcmp(option, 'parallel')
                    parfor i = 1:n
                        S(i) = obj.reach_star_single_input(in_sets(i), A);
                    end
                else
                    for i = 1:n
                        S(i) = obj.reach_star_single_input(in_sets(i), A);
                    end
                end
            end
        end

    end


    methods % helper methods

        function obj = setWeights(obj, W_node, W_edge, b)
            % Set weights with validation

            if size(W_node, 1) ~= size(W_edge, 1) || size(W_node, 2) ~= size(W_edge, 2)
                error('NodeWeights and EdgeWeights must have same dimensions');
            end
            if size(W_node, 2) ~= size(b, 1)
                error('Inconsistent dimensions: size(W,2) must equal size(b,1)');
            end
            if size(b, 2) ~= 1
                error('Bias vector should have one column');
            end
            obj.InputSize = size(W_node, 1);
            obj.OutputSize = size(W_node, 2);
            obj.NodeWeights = W_node;
            obj.EdgeWeights = W_edge;
            obj.Bias = b;
        end


        function obj = toGPU(obj)
            obj.NodeWeights = gpuArray(obj.NodeWeights);
            obj.EdgeWeights = gpuArray(obj.EdgeWeights);
            if ~isempty(obj.Bias)
                obj.Bias = gpuArray(obj.Bias);
            end
        end

        function obj = toCPU(obj)
            obj.NodeWeights = gather(obj.NodeWeights);
            obj.EdgeWeights = gather(obj.EdgeWeights);
            if ~isempty(obj.Bias)
                obj.Bias = gather(obj.Bias);
            end
        end

        function obj = changeParamsPrecision(obj, precision)
            if strcmp(precision, 'double')
                obj.NodeWeights = double(obj.NodeWeights);
                obj.EdgeWeights = double(obj.EdgeWeights);
                if ~isempty(obj.Bias)
                    obj.Bias = double(obj.Bias);
                end
            elseif strcmp(precision, 'single')
                obj.NodeWeights = single(obj.NodeWeights);
                obj.EdgeWeights = single(obj.EdgeWeights);
                if ~isempty(obj.Bias)
                    obj.Bias = single(obj.Bias);
                end
            else
                error('Precision must be "single" or "double"');
            end
        end

    end

end
