classdef HuGINEConvLayer < handle
    % HuGINEConvLayer - GIN with edge features from Hu et al. ICLR 2020
    %   "Strategies for Pre-training Graph Neural Networks"
    %
    %   Message passing (NO message ReLU):
    %     e_proj = E * W_edge + b_edge           (project edges to F_in)
    %     msg_ij = x_j + e_proj_ij               (linear message, no ReLU)
    %     agg_i = Σ_j w_ij * msg_ij              (weighted aggregation)
    %
    %   Self-loops are handled structurally: edge list includes (v,v) edges
    %   with zero edge features, so agg already includes the self term.
    %   No (1+eps) scaling.
    %
    %   MLP (applied to aggregation):
    %     h = agg * W1 + b1                      (MLP layer 1)
    %     h = ReLU(h)                             (MLP activation)
    %     Y = h * W2 + b2                         (MLP layer 2)
    %
    %   The outer ReLU (for intermediate layers) is handled by a separate
    %   ReluLayer in the GNN layer list, not inside this class.
    %
    %   Where:
    %     X: Input node features (N x F_in)
    %     E: Edge features (m x E_in), including self-loop edges
    %     adj_list: Edge list (m x 2) with [src, dst] pairs, including self-loops
    %     W_edge: Edge projection weights (E_in x F_in)
    %     W1: MLP first layer weights (F_in x hidden_dim)
    %     W2: MLP second layer weights (hidden_dim x F_out)
    %     w_ij: Edge weights for aggregation (optional)
    %     Y: Output node features (N x F_out)
    %
    % Main references:
    % 1) Hu et al., "Strategies for Pre-training Graph Neural Networks",
    %    ICLR 2020. https://arxiv.org/abs/1905.12265
    % 2) Xu et al., "How Powerful are Graph Neural Networks?",
    %    ICLR 2019. https://arxiv.org/abs/1810.00826
    %
    % Author: Anne Tumlin
    % Date: 03/17/2026

    properties
        Name = 'hu_gine_conv_layer';
        InputSize = 0;       % F_in: number of input features per node
        OutputSize = 0;      % F_out: number of output features per node
        HiddenSize = 0;      % hidden_dim: MLP hidden layer size
        EdgeInputSize = 0;   % E_in: number of input features per edge

        % MLP weights (2-layer MLP applied after message aggregation)
        MLPWeights1 = [];    % W1: F_in x hidden_dim
        MLPBias1 = [];       % b1: hidden_dim x 1
        MLPWeights2 = [];    % W2: hidden_dim x F_out
        MLPBias2 = [];       % b2: F_out x 1

        % Edge projection weights
        EdgeProjWeights = [];  % W_edge: E_in x F_in (projects edge to input dim)
        EdgeProjBias = [];     % b_edge: F_in x 1

        % Standard layer interface
        NumInputs = 1;
        InputNames = {'in'};
        NumOutputs = 1;
        OutputNames = {'out'};

        % Reach method (set by reach(), used by internal ReLU)
        reachMethod = 'approx-star';
    end


    methods % main methods

        % constructor of the class
        function obj = HuGINEConvLayer(varargin)
            % HuGINEConvLayer constructor
            % Usage:
            %   HuGINEConvLayer() - empty layer
            %   HuGINEConvLayer(W1, b1, W2, b2, W_edge, b_edge) - with weights
            %   HuGINEConvLayer(name, W1, b1, W2, b2, W_edge, b_edge) - with name
            %
            % @W1: MLP first layer weights (F_in x hidden_dim)
            % @b1: MLP first layer bias (hidden_dim x 1)
            % @W2: MLP second layer weights (hidden_dim x F_out)
            % @b2: MLP second layer bias (F_out x 1)
            % @W_edge: Edge projection weights (E_in x F_in)
            % @b_edge: Edge projection bias (F_in x 1)

            switch nargin

                case 7
                    name = varargin{1};
                    W1 = varargin{2};
                    b1 = varargin{3};
                    W2 = varargin{4};
                    b2 = varargin{5};
                    W_edge = varargin{6};
                    b_edge = varargin{7};

                    obj = obj.initWithWeights(name, W1, b1, W2, b2, W_edge, b_edge);

                case 6
                    W1 = varargin{1};
                    b1 = varargin{2};
                    W2 = varargin{3};
                    b2 = varargin{4};
                    W_edge = varargin{5};
                    b_edge = varargin{6};

                    obj = obj.initWithWeights('hu_gine_conv_layer', W1, b1, W2, b2, W_edge, b_edge);

                case 0
                    obj.Name = 'hu_gine_conv_layer';

                otherwise
                    error('Invalid number of inputs (should be 0, 6, or 7)');
            end

        end

        % helper for constructor
        function obj = initWithWeights(obj, name, W1, b1, W2, b2, W_edge, b_edge)
            if ~ischar(name) && ~isstring(name)
                error('Name must be a string or char array');
            end
            obj.Name = name;

            % Validate MLP layer 1: W1 is F_in x hidden_dim
            if size(W1, 2) ~= size(b1, 1)
                error('Inconsistent MLP layer 1 dimensions: size(W1,2) must equal size(b1,1)');
            end
            if size(b1, 2) ~= 1
                error('MLP bias 1 should have one column');
            end

            % Validate MLP layer 2: W2 is hidden_dim x F_out
            if size(W2, 1) ~= size(W1, 2)
                error('MLP layer 2 input must match layer 1 output: size(W2,1) must equal size(W1,2)');
            end
            if size(W2, 2) ~= size(b2, 1)
                error('Inconsistent MLP layer 2 dimensions: size(W2,2) must equal size(b2,1)');
            end
            if size(b2, 2) ~= 1
                error('MLP bias 2 should have one column');
            end

            % Validate edge weights: W_edge projects E_in to F_in
            if size(W_edge, 2) ~= size(W1, 1)
                error('Edge output dimension must match input dimension: size(W_edge,2) must equal size(W1,1)');
            end
            if size(b_edge, 1) ~= size(W_edge, 2)
                error('Edge bias dimension must match edge output: size(b_edge,1) must equal size(W_edge,2)');
            end
            if size(b_edge, 2) ~= 1
                error('Edge bias should have one column');
            end

            obj.InputSize = size(W1, 1);
            obj.HiddenSize = size(W1, 2);
            obj.OutputSize = size(W2, 2);
            obj.EdgeInputSize = size(W_edge, 1);
            obj.MLPWeights1 = W1;
            obj.MLPBias1 = b1;
            obj.MLPWeights2 = W2;
            obj.MLPBias2 = b2;
            obj.EdgeProjWeights = W_edge;
            obj.EdgeProjBias = b_edge;
        end

        % evaluation method
        function Y = evaluate(obj, X, E, adj_list, edge_weights)
            % Forward pass through Hu et al. GIN+E layer
            % @X: Input node features (N x F_in)
            % @E: Edge features (m x E_in), includes self-loop edges
            % @adj_list: Edge list (m x 2) with [src, dst] pairs, includes self-loops
            % @edge_weights: Normalized adjacency weights per edge (m x 1)
            %                Optional - defaults to ones if not provided
            % @Y: Output node features (N x F_out)
            %
            % Architecture (Hu et al. ICLR 2020):
            %   1. Project edge features to F_in
            %   2. Gather source node features to edges
            %   3. Add edge projection (NO ReLU)
            %   4. Aggregate to destination nodes (weighted, includes self-loops)
            %   5. MLP: Linear -> ReLU -> Linear

            if size(X, 2) ~= obj.InputSize
                error('Input feature dimension mismatch: expected %d, got %d', obj.InputSize, size(X, 2));
            end

            if size(E, 2) ~= obj.EdgeInputSize
                error('Edge feature dimension mismatch: expected %d, got %d', obj.EdgeInputSize, size(E, 2));
            end

            numNodes = size(X, 1);
            numEdges = size(adj_list, 1);

            if size(E, 1) ~= numEdges
                error('Number of edge features must match number of edges');
            end

            % Default edge weights to 1 if not provided
            if nargin < 5 || isempty(edge_weights)
                edge_weights = ones(numEdges, 1);
            end
            edge_weights = edge_weights(:);

            src_nodes = adj_list(:, 1);
            dst_nodes = adj_list(:, 2);

            % (1) Project edge features to F_in
            E_proj = E * obj.EdgeProjWeights;  % [m x F_in]
            if ~isempty(obj.EdgeProjBias)
                E_proj = E_proj + obj.EdgeProjBias';
            end

            % (2) Gather source node features to edges
            X_src = X(src_nodes, :);  % [m x F_in]

            % (3) Combine: NO ReLU (key difference from GINEConvLayer)
            edge_msg = X_src + E_proj;  % [m x F_in]

            % (4) Aggregate messages to destination nodes (weighted)
            %     Self-loops are in adj_list, so agg includes x_v + bias
            agg = zeros(numNodes, obj.InputSize);
            for e = 1:numEdges
                agg(dst_nodes(e), :) = agg(dst_nodes(e), :) + edge_weights(e) * edge_msg(e, :);
            end

            % (5) MLP: Linear -> ReLU -> Linear
            H = agg * obj.MLPWeights1;  % [N x hidden_dim]
            if ~isempty(obj.MLPBias1)
                H = H + obj.MLPBias1';
            end
            H = max(0, H);  % ReLU
            Y = H * obj.MLPWeights2;  % [N x F_out]
            if ~isempty(obj.MLPBias2)
                Y = Y + obj.MLPBias2';
            end

        end

        % main reachability analysis function
        function S = reach(varargin)
            % Reachability analysis for HuGINEConvLayer
            % @in_set: Input GraphStar set
            % @E: Edge features (matrix for node-only, Star for edge perturbation)
            % @adj_list: Edge list (m x 2), includes self-loops
            % @method: 'approx-star', 'exact-star', 'abs-dom'
            % @option: 'single' or 'parallel'
            % @edge_weights: Edge weights for aggregation (m x 1) - optional

            edge_weights = [];

            switch nargin

                case 9
                    obj = varargin{1};
                    in_sets = varargin{2};
                    E = varargin{3};
                    adj_list = varargin{4};
                    method = varargin{5};
                    option = varargin{6};
                    % relaxFactor = varargin{7};
                    % lp_solver = varargin{8};
                    edge_weights = varargin{9};

                case 8
                    obj = varargin{1};
                    in_sets = varargin{2};
                    E = varargin{3};
                    adj_list = varargin{4};
                    method = varargin{5};
                    option = varargin{6};

                case 7
                    obj = varargin{1};
                    in_sets = varargin{2};
                    E = varargin{3};
                    adj_list = varargin{4};
                    method = varargin{5};
                    option = varargin{6};

                case 6
                    obj = varargin{1};
                    in_sets = varargin{2};
                    E = varargin{3};
                    adj_list = varargin{4};
                    method = varargin{5};
                    option = varargin{6};

                case 5
                    obj = varargin{1};
                    in_sets = varargin{2};
                    E = varargin{3};
                    adj_list = varargin{4};
                    method = varargin{5};
                    option = [];

                case 4
                    obj = varargin{1};
                    in_sets = varargin{2};
                    E = varargin{3};
                    adj_list = varargin{4};
                    method = 'approx-star';
                    option = [];

                otherwise
                    error('Invalid number of input arguments (should be 3-8)');
            end

            obj.reachMethod = method;
            if strcmp(method, 'approx-star') || strcmp(method, 'exact-star') || strcmp(method, 'abs-dom') || contains(method, "relax-star")
                S = obj.reach_star_multipleInputs(in_sets, E, adj_list, option, edge_weights);
            elseif strcmp(method, 'approx-zono')
                error('Zonotope reachability for HuGINEConvLayer not yet implemented');
            else
                error('Unknown reachability method: %s', method);
            end

        end

    end


    methods % reachability methods

        % reachability for a single GraphStar input
        function gs_out = reach_star_single_input(obj, in_gs, E, adj_list, edge_weights)
            % Reachability through Hu et al. GIN+E layer for single GraphStar
            % @in_gs: Input GraphStar set
            % @E: Edge features (matrix for node-only mode)
            % @adj_list: Edge list (m x 2), includes self-loops
            % @edge_weights: Edge weights for aggregation (m x 1) - optional
            % @gs_out: Output GraphStar set

            if nargin < 5
                edge_weights = [];
            end

            % Type detection for perturbation mode
            if isa(E, 'Star') || isa(E, 'ImageStar') || isa(E, 'GraphStar')
                gs_out = obj.reach_with_edge_perturbation(in_gs, E, adj_list, edge_weights);
            else
                gs_out = obj.reach_node_only(in_gs, E, adj_list, edge_weights);
            end

        end

        % node-only perturbation reachability
        function gs_out = reach_node_only(obj, in_gs, E, adj_list, edge_weights)
            % Node-only reachability (edge features as constants)
            % Hu et al. GIN+E architecture:
            %   1. Project edge features to F_in (constant)
            %   2. Gather source features to edges
            %   3. Add edge projection (NO ReLU — purely linear)
            %   4. Aggregate to nodes (weighted, includes self-loops)
            %   5. MLP layer 1 (linear)
            %   6. MLP ReLU (only ReLU in the layer)
            %   7. MLP layer 2 (linear)
            %
            % @in_gs: Input GraphStar set
            % @E: Edge features (m x E_in) - constant matrix
            % @adj_list: Edge list (m x 2), includes self-loops
            % @edge_weights: Edge weights for aggregation (m x 1) - optional
            % @gs_out: Output GraphStar set

            if ~isa(in_gs, 'GraphStar')
                error('Input must be a GraphStar');
            end

            if in_gs.numFeatures ~= obj.InputSize
                error('Input feature dimension mismatch: expected %d, got %d', obj.InputSize, in_gs.numFeatures);
            end

            numNodes = in_gs.numNodes;
            numEdges = size(adj_list, 1);

            if nargin < 5 || isempty(edge_weights)
                edge_weights = ones(numEdges, 1);
            end
            edge_weights = edge_weights(:);

            src_nodes = adj_list(:, 1);
            dst_nodes = adj_list(:, 2);

            % --- Message passing phase (entirely linear, no ReLU) ---

            % (1) Project edge features to F_in (constant, added to center only)
            E_proj = E * obj.EdgeProjWeights;  % [m x F_in]
            if ~isempty(obj.EdgeProjBias)
                E_proj = E_proj + obj.EdgeProjBias';
            end

            % (2) Gather source node features to edges
            V_edge = obj.gather_src_features_star(in_gs.V, src_nodes);  % [m x F_in x K]

            % (3) Add constant edge projection to center (NO ReLU)
            V_edge(:, :, 1) = V_edge(:, :, 1) + E_proj;

            % (4) Aggregate edge messages back to nodes (weighted)
            %     Self-loops already in adj_list, so this includes self-connection
            V_agg = obj.aggregate_edges_to_nodes_weighted_star(V_edge, dst_nodes, numNodes, edge_weights);

            % Constraints unchanged through the linear message phase
            C_agg = in_gs.C;
            d_agg = in_gs.d;
            pred_lb_agg = in_gs.pred_lb;
            pred_ub_agg = in_gs.pred_ub;

            % --- MLP phase ---

            % (5) MLP layer 1: linear transform V_agg * W1 + b1
            V_mlp = obj.map_features_star(V_agg, obj.MLPWeights1);  % [N x hidden x K]
            if ~isempty(obj.MLPBias1)
                V_mlp(:, :, 1) = V_mlp(:, :, 1) + obj.MLPBias1';
            end

            % (6) MLP ReLU: flatten to Star, apply ReLU, reshape back
            [N, H_dim, K_mlp] = size(V_mlp);
            V_flat_mlp = reshape(permute(V_mlp, [2, 1, 3]), [N * H_dim, K_mlp]);

            mlp_star = Star(V_flat_mlp, C_agg, d_agg, pred_lb_agg, pred_ub_agg);

            L_relu = ReluLayer();
            mlp_star_out = L_relu.reach(mlp_star, obj.reachMethod);
            if iscell(mlp_star_out)
                mlp_star_out = mlp_star_out{1};
            end

            K_mlp_out = mlp_star_out.nVar + 1;
            V_flat_mlp_out = mlp_star_out.V;
            V_mlp_relu = permute(reshape(V_flat_mlp_out, [H_dim, N, K_mlp_out]), [2, 1, 3]);

            C_out = mlp_star_out.C;
            d_out = mlp_star_out.d;
            pred_lb_out = mlp_star_out.predicate_lb;
            pred_ub_out = mlp_star_out.predicate_ub;

            % (7) MLP layer 2: linear transform
            V_out = obj.map_features_star(V_mlp_relu, obj.MLPWeights2);  % [N x F_out x K]
            if ~isempty(obj.MLPBias2)
                V_out(:, :, 1) = V_out(:, :, 1) + obj.MLPBias2';
            end

            % Create output GraphStar
            gs_out = GraphStar(V_out, C_out, d_out, pred_lb_out, pred_ub_out);

        end

        % edge perturbation reachability
        function gs_out = reach_with_edge_perturbation(obj, in_gs, E_star, adj_list, edge_weights)
            % Edge perturbation reachability (both node and edge features perturbed)
            % @in_gs: Input GraphStar set for node features
            % @E_star: Edge features as Star, ImageStar, or GraphStar
            % @adj_list: Edge list (m x 2), includes self-loops
            % @edge_weights: Edge weights for aggregation (m x 1) - optional
            % @gs_out: Output GraphStar set

            if ~isa(in_gs, 'GraphStar')
                error('Input must be a GraphStar');
            end

            if in_gs.numFeatures ~= obj.InputSize
                error('Input feature dimension mismatch: expected %d, got %d', obj.InputSize, in_gs.numFeatures);
            end

            numNodes = in_gs.numNodes;
            numEdges = size(adj_list, 1);

            if nargin < 5 || isempty(edge_weights)
                edge_weights = ones(numEdges, 1);
            end
            edge_weights = edge_weights(:);

            src_nodes = adj_list(:, 1);
            dst_nodes = adj_list(:, 2);

            % Extract edge Star properties
            if isa(E_star, 'ImageStar')
                E_V = squeeze(E_star.V);
                if ndims(E_V) == 2
                    E_V = reshape(E_V, [size(E_V, 1), size(E_V, 2), 1]);
                end
                E_C = E_star.C;
                E_d = E_star.d;
                E_pred_lb = E_star.pred_lb;
                E_pred_ub = E_star.pred_ub;
            elseif isa(E_star, 'GraphStar')
                E_V = E_star.V;
                E_C = E_star.C;
                E_d = E_star.d;
                E_pred_lb = E_star.pred_lb;
                E_pred_ub = E_star.pred_ub;
            else
                E_V = reshape(E_star.V, [numEdges, obj.EdgeInputSize, E_star.nVar + 1]);
                E_C = E_star.C;
                E_d = E_star.d;
                E_pred_lb = E_star.predicate_lb;
                E_pred_ub = E_star.predicate_ub;
            end

            % (1) Project edge features to F_in (edge features are uncertain)
            E_proj_V = obj.map_features_star(E_V, obj.EdgeProjWeights);  % [m x F_in x K_edge]
            if ~isempty(obj.EdgeProjBias)
                E_proj_V(:, :, 1) = E_proj_V(:, :, 1) + obj.EdgeProjBias';
            end

            % (2) Gather source node features to edges
            V_node_edge = obj.gather_src_features_star(in_gs.V, src_nodes);  % [m x F_in x K_node]

            % (3) Combine via Minkowski sum (direct, no padding)
            C_node = in_gs.C; d_node = in_gs.d;
            lb_node = in_gs.pred_lb; ub_node = in_gs.pred_ub;
            C_edge = E_C; d_edge = E_d;
            lb_edge = E_pred_lb; ub_edge = E_pred_ub;

            K_node = size(V_node_edge, 3);
            K_edge = size(E_proj_V, 3);

            % Normalize empty constraints for blkdiag
            if isempty(C_node), C_node = zeros(0, max(K_node-1, 0)); end
            if isempty(C_edge), C_edge = zeros(0, max(K_edge-1, 0)); end
            if isempty(d_node), d_node = zeros(0, 1); end
            if isempty(d_edge), d_edge = zeros(0, 1); end
            if isempty(lb_node), lb_node = -ones(max(K_node-1, 0), 1); end
            if isempty(ub_node), ub_node = ones(max(K_node-1, 0), 1); end
            if isempty(lb_edge), lb_edge = -ones(max(K_edge-1, 0), 1); end
            if isempty(ub_edge), ub_edge = ones(max(K_edge-1, 0), 1); end

            % NO ReLU here — messages are purely linear
            % Combined message: centers add, generators concatenate
            V_combined_msg = zeros(numEdges, obj.InputSize, K_node + K_edge - 1, 'like', V_node_edge);
            V_combined_msg(:, :, 1) = V_node_edge(:, :, 1) + E_proj_V(:, :, 1);
            V_combined_msg(:, :, 2:K_node) = V_node_edge(:, :, 2:end);
            V_combined_msg(:, :, K_node+1:end) = E_proj_V(:, :, 2:end);

            C_combined = blkdiag(C_node, C_edge);
            d_combined = [d_node; d_edge];
            pred_lb_combined = [lb_node; lb_edge];
            pred_ub_combined = [ub_node; ub_edge];

            % (4) Aggregate edge messages back to nodes (weighted)
            K_msg = K_node + K_edge - 1;
            V_agg = obj.aggregate_edges_to_nodes_weighted_star(V_combined_msg, dst_nodes, numNodes, edge_weights);

            % --- MLP phase ---

            % (5) MLP layer 1
            V_mlp = obj.map_features_star(V_agg, obj.MLPWeights1);
            if ~isempty(obj.MLPBias1)
                V_mlp(:, :, 1) = V_mlp(:, :, 1) + obj.MLPBias1';
            end

            % (6) MLP ReLU
            [N, H_dim, K_mlp] = size(V_mlp);
            V_flat_mlp = reshape(permute(V_mlp, [2, 1, 3]), [N * H_dim, K_mlp]);

            mlp_star = Star(V_flat_mlp, C_combined, d_combined, pred_lb_combined, pred_ub_combined);

            L_relu = ReluLayer();
            mlp_star_out = L_relu.reach(mlp_star, obj.reachMethod);
            if iscell(mlp_star_out)
                mlp_star_out = mlp_star_out{1};
            end

            K_mlp_out = mlp_star_out.nVar + 1;
            V_flat_mlp_out = mlp_star_out.V;
            V_mlp_relu = permute(reshape(V_flat_mlp_out, [H_dim, N, K_mlp_out]), [2, 1, 3]);

            C_out = mlp_star_out.C;
            d_out = mlp_star_out.d;
            pred_lb_out = mlp_star_out.predicate_lb;
            pred_ub_out = mlp_star_out.predicate_ub;

            % (7) MLP layer 2
            V_out = obj.map_features_star(V_mlp_relu, obj.MLPWeights2);
            if ~isempty(obj.MLPBias2)
                V_out(:, :, 1) = V_out(:, :, 1) + obj.MLPBias2';
            end

            gs_out = GraphStar(V_out, C_out, d_out, pred_lb_out, pred_ub_out);

        end

        % reachability for multiple GraphStar inputs
        function S = reach_star_multipleInputs(obj, in_sets, E, adj_list, option, edge_weights)

            if nargin < 6
                edge_weights = [];
            end

            n = length(in_sets);

            if n == 1
                S = obj.reach_star_single_input(in_sets, E, adj_list, edge_weights);
            else
                S(n) = GraphStar;

                if strcmp(option, 'parallel')
                    parfor i = 1:n
                        S(i) = obj.reach_star_single_input(in_sets(i), E, adj_list, edge_weights);
                    end
                else
                    for i = 1:n
                        S(i) = obj.reach_star_single_input(in_sets(i), E, adj_list, edge_weights);
                    end
                end
            end

        end

    end


    methods (Access = private) % private helper methods

        % Gather source node features to edge representation
        function V_edge = gather_src_features_star(~, V_node, src_nodes)
            % V_node: [N x F x K] -> V_edge: [m x F x K]
            [~, F, K] = size(V_node);
            m = numel(src_nodes);
            V_edge = zeros(m, F, K, 'like', V_node);
            for k = 1:K
                V_edge(:, :, k) = V_node(src_nodes, :, k);
            end
        end

        % Aggregate edge features back to nodes with edge weights
        function V_node = aggregate_edges_to_nodes_weighted_star(~, V_edge, dst_nodes, numNodes, edge_weights)
            % V_edge: [m x F x K] -> V_node: [N x F x K]
            [m, F, K] = size(V_edge);
            V_node = zeros(numNodes, F, K, 'like', V_edge);
            for k = 1:K
                for e = 1:m
                    V_node(dst_nodes(e), :, k) = V_node(dst_nodes(e), :, k) + edge_weights(e) * V_edge(e, :, k);
                end
            end
        end

        % Apply linear transformation to all generators
        function V_out = map_features_star(~, V_in, W)
            % V_in: [N x F_in x K], W: [F_in x F_out]
            % V_out: [N x F_out x K]
            [N, ~, K] = size(V_in);
            F_out = size(W, 2);
            V_out = zeros(N, F_out, K, 'like', V_in);
            for k = 1:K
                V_out(:, :, k) = V_in(:, :, k) * W;
            end
        end

    end


    methods % helper methods

        function obj = toGPU(obj)
            obj.MLPWeights1 = gpuArray(obj.MLPWeights1);
            obj.MLPWeights2 = gpuArray(obj.MLPWeights2);
            obj.EdgeProjWeights = gpuArray(obj.EdgeProjWeights);
            if ~isempty(obj.MLPBias1), obj.MLPBias1 = gpuArray(obj.MLPBias1); end
            if ~isempty(obj.MLPBias2), obj.MLPBias2 = gpuArray(obj.MLPBias2); end
            if ~isempty(obj.EdgeProjBias), obj.EdgeProjBias = gpuArray(obj.EdgeProjBias); end
        end

        function obj = toCPU(obj)
            obj.MLPWeights1 = gather(obj.MLPWeights1);
            obj.MLPWeights2 = gather(obj.MLPWeights2);
            obj.EdgeProjWeights = gather(obj.EdgeProjWeights);
            if ~isempty(obj.MLPBias1), obj.MLPBias1 = gather(obj.MLPBias1); end
            if ~isempty(obj.MLPBias2), obj.MLPBias2 = gather(obj.MLPBias2); end
            if ~isempty(obj.EdgeProjBias), obj.EdgeProjBias = gather(obj.EdgeProjBias); end
        end

        function obj = changeParamsPrecision(obj, precision)
            if strcmp(precision, 'double')
                obj.MLPWeights1 = double(obj.MLPWeights1);
                obj.MLPWeights2 = double(obj.MLPWeights2);
                obj.EdgeProjWeights = double(obj.EdgeProjWeights);
                if ~isempty(obj.MLPBias1), obj.MLPBias1 = double(obj.MLPBias1); end
                if ~isempty(obj.MLPBias2), obj.MLPBias2 = double(obj.MLPBias2); end
                if ~isempty(obj.EdgeProjBias), obj.EdgeProjBias = double(obj.EdgeProjBias); end
            elseif strcmp(precision, 'single')
                obj.MLPWeights1 = single(obj.MLPWeights1);
                obj.MLPWeights2 = single(obj.MLPWeights2);
                obj.EdgeProjWeights = single(obj.EdgeProjWeights);
                if ~isempty(obj.MLPBias1), obj.MLPBias1 = single(obj.MLPBias1); end
                if ~isempty(obj.MLPBias2), obj.MLPBias2 = single(obj.MLPBias2); end
                if ~isempty(obj.EdgeProjBias), obj.EdgeProjBias = single(obj.EdgeProjBias); end
            else
                error('Precision must be "single" or "double"');
            end
        end

    end

end
