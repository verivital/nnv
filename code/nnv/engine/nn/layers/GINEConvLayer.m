classdef GINEConvLayer < handle
    % The GINEConvLayer class for full Graph Isomorphism Networks with Edge features
    %   Matches PyTorch Geometric's GINEConv architecture with 2-layer MLPs.
    %
    %   Message passing:
    %     e_proj = E * W_edge + b_edge           (project edges to F_in)
    %     msg_ij = ReLU(x_j + e_proj_ij)         (edge message with ReLU)
    %     agg_i = Σ_j w_ij * msg_ij              (weighted aggregation)
    %     combined_i = (1+eps)*x_i + agg_i       (self-loop + aggregation)
    %
    %   MLP (applied to combined):
    %     h = combined * W1 + b1                  (MLP layer 1)
    %     h = ReLU(h)                             (MLP activation)
    %     Y = h * W2 + b2                         (MLP layer 2)
    %
    %   Where:
    %     X: Input node features (N x F_in)
    %     E: Edge features (m x E_in)
    %     adj_list: Edge list (m x 2) with [src, dst] pairs
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
    % Date: 02/11/2026

    properties
        Name = 'gine_conv_layer';
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
        EdgeWeights = [];    % W_edge: E_in x F_in (projects edge to input dim)
        EdgeBias = [];       % b_edge: F_in x 1

        % Epsilon parameter for self-loop scaling
        Epsilon = 0;

        % Standard layer interface
        NumInputs = 1;
        InputNames = {'in'};
        NumOutputs = 1;
        OutputNames = {'out'};
    end


    methods % main methods

        % constructor of the class
        function obj = GINEConvLayer(varargin)
            % GINEConvLayer constructor
            % Usage:
            %   GINEConvLayer() - empty layer
            %   GINEConvLayer(W1, b1, W2, b2, W_edge, b_edge) - with weights
            %   GINEConvLayer(name, W1, b1, W2, b2, W_edge, b_edge) - with name
            %   GINEConvLayer(name, W1, b1, W2, b2, W_edge, b_edge, epsilon)
            %
            % @W1: MLP first layer weights (F_in x hidden_dim)
            % @b1: MLP first layer bias (hidden_dim x 1)
            % @W2: MLP second layer weights (hidden_dim x F_out)
            % @b2: MLP second layer bias (F_out x 1)
            % @W_edge: Edge projection weights (E_in x F_in)
            % @b_edge: Edge projection bias (F_in x 1)
            % @epsilon: Self-loop scaling factor

            switch nargin

                case 8
                    name = varargin{1};
                    W1 = varargin{2};
                    b1 = varargin{3};
                    W2 = varargin{4};
                    b2 = varargin{5};
                    W_edge = varargin{6};
                    b_edge = varargin{7};
                    epsilon = varargin{8};

                    obj = obj.initWithWeights(name, W1, b1, W2, b2, W_edge, b_edge, epsilon);

                case 7
                    name = varargin{1};
                    W1 = varargin{2};
                    b1 = varargin{3};
                    W2 = varargin{4};
                    b2 = varargin{5};
                    W_edge = varargin{6};
                    b_edge = varargin{7};

                    obj = obj.initWithWeights(name, W1, b1, W2, b2, W_edge, b_edge, 0);

                case 6
                    W1 = varargin{1};
                    b1 = varargin{2};
                    W2 = varargin{3};
                    b2 = varargin{4};
                    W_edge = varargin{5};
                    b_edge = varargin{6};

                    obj = obj.initWithWeights('gine_conv_layer', W1, b1, W2, b2, W_edge, b_edge, 0);

                case 0
                    obj.Name = 'gine_conv_layer';

                otherwise
                    error('Invalid number of inputs (should be 0, 6, 7, or 8)');
            end

        end

        % helper for constructor
        function obj = initWithWeights(obj, name, W1, b1, W2, b2, W_edge, b_edge, epsilon)
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
            obj.EdgeWeights = W_edge;
            obj.EdgeBias = b_edge;
            obj.Epsilon = epsilon;
        end

        % evaluation method
        function Y = evaluate(obj, X, E, adj_list, edge_weights)
            % Forward pass through GINEConv layer (full GINE with MLP)
            % @X: Input node features (N x F_in)
            % @E: Edge features (m x E_in)
            % @adj_list: Edge list (m x 2) with [src, dst] pairs
            % @edge_weights: Normalized adjacency weights per edge (m x 1)
            %                Optional - defaults to ones if not provided
            % @Y: Output node features (N x F_out)
            %
            % Architecture (matches PyTorch Geometric GINEConv):
            %   1. Project edge features to F_in
            %   2. Gather source node features to edges
            %   3. Add edge projection and apply ReLU
            %   4. Aggregate to destination nodes (weighted)
            %   5. Self-loop: (1+eps)*x + agg
            %   6. MLP: Linear -> ReLU -> Linear

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
            E_proj = E * obj.EdgeWeights;  % [m x F_in]
            if ~isempty(obj.EdgeBias)
                E_proj = E_proj + obj.EdgeBias';
            end

            % (2) Gather source node features to edges
            X_src = X(src_nodes, :);  % [m x F_in]

            % (3) Combine and ReLU
            edge_msg = max(0, X_src + E_proj);  % [m x F_in]

            % (4) Aggregate messages to destination nodes (weighted)
            agg = zeros(numNodes, obj.InputSize);
            for e = 1:numEdges
                agg(dst_nodes(e), :) = agg(dst_nodes(e), :) + edge_weights(e) * edge_msg(e, :);
            end

            % (5) Self-loop: (1+eps)*x + agg
            combined = (1 + obj.Epsilon) * X + agg;  % [N x F_in]

            % (6) MLP: Linear -> ReLU -> Linear
            H = combined * obj.MLPWeights1;  % [N x hidden_dim]
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
            % Reachability analysis for GINEConvLayer
            % @in_set: Input GraphStar set
            % @E: Edge features (matrix for node-only, Star for edge perturbation)
            % @adj_list: Edge list (m x 2)
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

            if strcmp(method, 'approx-star') || strcmp(method, 'exact-star') || strcmp(method, 'abs-dom') || contains(method, "relax-star")
                S = obj.reach_star_multipleInputs(in_sets, E, adj_list, option, edge_weights);
            elseif strcmp(method, 'approx-zono')
                error('Zonotope reachability for GINEConvLayer not yet implemented');
            else
                error('Unknown reachability method: %s', method);
            end

        end

    end


    methods % reachability methods

        % reachability for a single GraphStar input
        function gs_out = reach_star_single_input(obj, in_gs, E, adj_list, edge_weights)
            % Reachability through GINEConv layer for single GraphStar
            % @in_gs: Input GraphStar set
            % @E: Edge features (matrix for node-only mode)
            % @adj_list: Edge list (m x 2)
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
            % Full GINEConv architecture:
            %   1. Project edge features to F_in (constant)
            %   2. Gather source features to edges
            %   3. Add edge projection and ReLU (message ReLU)
            %   4. Aggregate to nodes (weighted)
            %   5. Self-loop: (1+eps)*x + agg
            %   6. MLP layer 1 (linear)
            %   7. MLP ReLU
            %   8. MLP layer 2 (linear)
            %
            % @in_gs: Input GraphStar set
            % @E: Edge features (m x E_in) - constant matrix
            % @adj_list: Edge list (m x 2)
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

            % --- Message passing phase ---

            % (1) Project edge features to F_in (constant, added to center only)
            E_proj = E * obj.EdgeWeights;  % [m x F_in]
            if ~isempty(obj.EdgeBias)
                E_proj = E_proj + obj.EdgeBias';
            end

            % (2) Gather source node features to edges
            V_edge = obj.gather_src_features_star(in_gs.V, src_nodes);  % [m x F_in x K]

            % (3) Add constant edge projection to center
            V_edge(:, :, 1) = V_edge(:, :, 1) + E_proj;

            % (4) Apply message ReLU using NNV's ReluLayer
            [m_edges, F_in, K_orig] = size(V_edge);
            V_flat = reshape(permute(V_edge, [2, 1, 3]), [m_edges * F_in, K_orig]);

            edge_star = Star(V_flat, in_gs.C, in_gs.d, in_gs.pred_lb, in_gs.pred_ub);

            L_relu = ReluLayer();
            edge_star_out = L_relu.reach(edge_star, 'approx-star');
            if iscell(edge_star_out)
                edge_star_out = edge_star_out{1};
            end

            K_msg = edge_star_out.nVar + 1;
            V_flat_out = edge_star_out.V;
            V_edge_relu = permute(reshape(V_flat_out, [F_in, m_edges, K_msg]), [2, 1, 3]);

            C_msg = edge_star_out.C;
            d_msg = edge_star_out.d;
            pred_lb_msg = edge_star_out.predicate_lb;
            pred_ub_msg = edge_star_out.predicate_ub;

            % (5) Aggregate edge messages back to nodes (weighted)
            V_agg = obj.aggregate_edges_to_nodes_weighted_star(V_edge_relu, dst_nodes, numNodes, edge_weights);

            % (6) Self-loop: (1+eps)*x + agg
            % Expand in_gs.V to match K_msg if needed
            K_in = size(in_gs.V, 3);
            if K_msg > K_in
                V_self = zeros(numNodes, F_in, K_msg, 'like', in_gs.V);
                V_self(:, :, 1:K_in) = in_gs.V;
            else
                V_self = in_gs.V;
            end
            V_combined = (1 + obj.Epsilon) * V_self + V_agg;

            % --- MLP phase ---
            % Convert to Star, apply MLP layers, convert back to GraphStar

            % MLP layer 1: linear transform V_combined * W1 + b1
            V_mlp = obj.map_features_star(V_combined, obj.MLPWeights1);  % [N x hidden x K]
            if ~isempty(obj.MLPBias1)
                V_mlp(:, :, 1) = V_mlp(:, :, 1) + obj.MLPBias1';
            end

            % MLP ReLU: flatten to Star, apply ReLU, reshape back
            [N, H_dim, K_mlp] = size(V_mlp);
            V_flat_mlp = reshape(permute(V_mlp, [2, 1, 3]), [N * H_dim, K_mlp]);

            mlp_star = Star(V_flat_mlp, C_msg, d_msg, pred_lb_msg, pred_ub_msg);

            mlp_star_out = L_relu.reach(mlp_star, 'approx-star');
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

            % MLP layer 2: linear transform
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
            % @adj_list: Edge list (m x 2)
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
            E_proj_V = obj.map_features_star(E_V, obj.EdgeWeights);  % [m x F_in x K_edge]
            if ~isempty(obj.EdgeBias)
                E_proj_V(:, :, 1) = E_proj_V(:, :, 1) + obj.EdgeBias';
            end

            % (2) Gather source node features to edges
            V_node_edge = obj.gather_src_features_star(in_gs.V, src_nodes);  % [m x F_in x K_node]

            % (3) Combine via Minkowski sum (direct, no padding)
            % Avoids inflating generator count from 2*K_pad-1 to K_node+K_edge-1
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

            V_combined_msg = zeros(numEdges, obj.InputSize, K_node + K_edge - 1, 'like', V_node_edge);
            V_combined_msg(:, :, 1) = V_node_edge(:, :, 1) + E_proj_V(:, :, 1);
            V_combined_msg(:, :, 2:K_node) = V_node_edge(:, :, 2:end);
            V_combined_msg(:, :, K_node+1:end) = E_proj_V(:, :, 2:end);

            C_combined = blkdiag(C_node, C_edge);
            d_combined = [d_node; d_edge];
            pred_lb_combined = [lb_node; lb_edge];
            pred_ub_combined = [ub_node; ub_edge];

            % (4) Apply message ReLU
            [m_edges, F_in, K_comb] = size(V_combined_msg);
            V_flat = reshape(permute(V_combined_msg, [2, 1, 3]), [m_edges * F_in, K_comb]);

            edge_star = Star(V_flat, C_combined, d_combined, pred_lb_combined, pred_ub_combined);

            L_relu = ReluLayer();
            edge_star_out = L_relu.reach(edge_star, 'approx-star');
            if iscell(edge_star_out)
                edge_star_out = edge_star_out{1};
            end

            K_msg = edge_star_out.nVar + 1;
            V_flat_out = edge_star_out.V;
            V_edge_relu = permute(reshape(V_flat_out, [F_in, m_edges, K_msg]), [2, 1, 3]);

            C_msg = edge_star_out.C;
            d_msg = edge_star_out.d;
            pred_lb_msg = edge_star_out.predicate_lb;
            pred_ub_msg = edge_star_out.predicate_ub;

            % (5) Aggregate edge messages back to nodes (weighted)
            V_agg = obj.aggregate_edges_to_nodes_weighted_star(V_edge_relu, dst_nodes, numNodes, edge_weights);

            % (6) Self-loop: (1+eps)*x + agg
            K_in = size(in_gs.V, 3);
            if K_msg > K_in
                V_self = zeros(numNodes, F_in, K_msg, 'like', in_gs.V);
                V_self(:, :, 1:K_in) = in_gs.V;
            else
                V_self = in_gs.V;
            end
            V_combined = (1 + obj.Epsilon) * V_self + V_agg;

            % --- MLP phase (same as node-only) ---

            % MLP layer 1
            V_mlp = obj.map_features_star(V_combined, obj.MLPWeights1);
            if ~isempty(obj.MLPBias1)
                V_mlp(:, :, 1) = V_mlp(:, :, 1) + obj.MLPBias1';
            end

            % MLP ReLU
            [N, H_dim, K_mlp] = size(V_mlp);
            V_flat_mlp = reshape(permute(V_mlp, [2, 1, 3]), [N * H_dim, K_mlp]);

            mlp_star = Star(V_flat_mlp, C_msg, d_msg, pred_lb_msg, pred_ub_msg);

            mlp_star_out = L_relu.reach(mlp_star, 'approx-star');
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

            % MLP layer 2
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

        % Pad two star representations to have same number of generators
        function [V1_out, C1_out, d1_out, lb1_out, ub1_out, ...
                  V2_out, C2_out, d2_out, lb2_out, ub2_out] = ...
                pad_to_same_K(~, V1, C1, d1, lb1, ub1, V2, C2, d2, lb2, ub2)

            K1 = size(V1, 3);
            K2 = size(V2, 3);

            if K1 == K2
                V1_out = V1; V2_out = V2;
                C1_out = C1; C2_out = C2;
                d1_out = d1; d2_out = d2;
                lb1_out = lb1; ub1_out = ub1;
                lb2_out = lb2; ub2_out = ub2;
                return;
            end

            if isempty(C1), C1 = zeros(0, K1 - 1); end
            if isempty(C2), C2 = zeros(0, K2 - 1); end
            if isempty(d1), d1 = zeros(0, 1); end
            if isempty(d2), d2 = zeros(0, 1); end
            if isempty(lb1), lb1 = -ones(K1 - 1, 1); end
            if isempty(ub1), ub1 = ones(K1 - 1, 1); end
            if isempty(lb2), lb2 = -ones(K2 - 1, 1); end
            if isempty(ub2), ub2 = ones(K2 - 1, 1); end

            K_max = max(K1, K2);

            if K1 < K_max
                V1_out = zeros(size(V1, 1), size(V1, 2), K_max, 'like', V1);
                V1_out(:, :, 1:K1) = V1;
                C1_out = [C1, zeros(size(C1, 1), K_max - K1)];
                d1_out = d1;
                lb1_out = [lb1; zeros(K_max - K1, 1)];
                ub1_out = [ub1; zeros(K_max - K1, 1)];
            else
                V1_out = V1; C1_out = C1; d1_out = d1;
                lb1_out = lb1; ub1_out = ub1;
            end

            if K2 < K_max
                V2_out = zeros(size(V2, 1), size(V2, 2), K_max, 'like', V2);
                V2_out(:, :, 1:K2) = V2;
                C2_out = [C2, zeros(size(C2, 1), K_max - K2)];
                d2_out = d2;
                lb2_out = [lb2; zeros(K_max - K2, 1)];
                ub2_out = [ub2; zeros(K_max - K2, 1)];
            else
                V2_out = V2; C2_out = C2; d2_out = d2;
                lb2_out = lb2; ub2_out = ub2;
            end

        end

    end


    methods % helper methods

        function obj = toGPU(obj)
            obj.MLPWeights1 = gpuArray(obj.MLPWeights1);
            obj.MLPWeights2 = gpuArray(obj.MLPWeights2);
            obj.EdgeWeights = gpuArray(obj.EdgeWeights);
            if ~isempty(obj.MLPBias1), obj.MLPBias1 = gpuArray(obj.MLPBias1); end
            if ~isempty(obj.MLPBias2), obj.MLPBias2 = gpuArray(obj.MLPBias2); end
            if ~isempty(obj.EdgeBias), obj.EdgeBias = gpuArray(obj.EdgeBias); end
        end

        function obj = toCPU(obj)
            obj.MLPWeights1 = gather(obj.MLPWeights1);
            obj.MLPWeights2 = gather(obj.MLPWeights2);
            obj.EdgeWeights = gather(obj.EdgeWeights);
            if ~isempty(obj.MLPBias1), obj.MLPBias1 = gather(obj.MLPBias1); end
            if ~isempty(obj.MLPBias2), obj.MLPBias2 = gather(obj.MLPBias2); end
            if ~isempty(obj.EdgeBias), obj.EdgeBias = gather(obj.EdgeBias); end
        end

        function obj = changeParamsPrecision(obj, precision)
            if strcmp(precision, 'double')
                obj.MLPWeights1 = double(obj.MLPWeights1);
                obj.MLPWeights2 = double(obj.MLPWeights2);
                obj.EdgeWeights = double(obj.EdgeWeights);
                if ~isempty(obj.MLPBias1), obj.MLPBias1 = double(obj.MLPBias1); end
                if ~isempty(obj.MLPBias2), obj.MLPBias2 = double(obj.MLPBias2); end
                if ~isempty(obj.EdgeBias), obj.EdgeBias = double(obj.EdgeBias); end
            elseif strcmp(precision, 'single')
                obj.MLPWeights1 = single(obj.MLPWeights1);
                obj.MLPWeights2 = single(obj.MLPWeights2);
                obj.EdgeWeights = single(obj.EdgeWeights);
                if ~isempty(obj.MLPBias1), obj.MLPBias1 = single(obj.MLPBias1); end
                if ~isempty(obj.MLPBias2), obj.MLPBias2 = single(obj.MLPBias2); end
                if ~isempty(obj.EdgeBias), obj.EdgeBias = single(obj.EdgeBias); end
            else
                error('Precision must be "single" or "double"');
            end
        end

    end

end
