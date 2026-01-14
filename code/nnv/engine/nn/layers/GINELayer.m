classdef GINELayer < handle
    % The GINELayer class for Graph Isomorphism Networks with Edge features
    %   Performs message passing with edge features:
    %   H = X * W_node                           (transform nodes if F_in != F_out)
    %   msg_ij = ReLU(h_i + E_ij * W_edge + b_edge)  (edge message)
    %   agg_j = Σ_i w_ij * msg_ij                (weighted aggregation)
    %   Y = H + agg + b_node                     (residual connection)
    %
    %   Where:
    %     X: Input node features (N x F_in)
    %     E: Edge features (m x E_in)
    %     adj_list: Edge list (m x 2) with [src, dst] pairs
    %     W_node: Node transformation weights (F_in x F_out)
    %     W_edge: Edge transformation weights (E_in x F_out)
    %     w_ij: Edge weights for aggregation (optional)
    %     Y: Output node features (N x F_out)
    %
    % Main references:
    % 1) Hu et al., "Strategies for Pre-training Graph Neural Networks",
    %    ICLR 2020. https://arxiv.org/abs/1905.12265
    % 2) GNNV: Graph Neural Network Verification prototype
    %    (internal reference implementation)
    % 3) NNV verification tool documentation
    %    https://github.com/verivital/nnv
    %
    % Author: Anne Tumlin
    % Date: 01/13/2026

    properties
        Name = 'gine_layer';
        InputSize = 0;       % F_in: number of input features per node
        OutputSize = 0;      % F_out: number of output features per node
        EdgeInputSize = 0;   % E_in: number of input features per edge

        % Node transformation weights
        Weights = [];        % W_node: F_in x F_out
        Bias = [];           % b_node: F_out x 1

        % Edge transformation weights
        EdgeWeights = [];    % W_edge: E_in x F_out (projects edge to output dim)
        EdgeBias = [];       % b_edge: F_out x 1

        % Epsilon parameter for self-loop
        Epsilon = 0;

        % Standard layer interface
        NumInputs = 1;
        InputNames = {'in'};
        NumOutputs = 1;
        OutputNames = {'out'};
    end


    methods % main methods

        % constructor of the class
        function obj = GINELayer(varargin)
            % GINELayer constructor
            % Usage:
            %   GINELayer() - empty layer
            %   GINELayer(W_node, b_node, W_edge, b_edge) - with weights
            %   GINELayer(name, W_node, b_node, W_edge, b_edge) - with name
            %   GINELayer(name, W_node, b_node, W_edge, b_edge, epsilon)
            %
            % @W_node: Node weight matrix (F_in x F_out)
            % @b_node: Node bias vector (F_out x 1)
            % @W_edge: Edge weight matrix (E_in x F_in)
            % @b_edge: Edge bias vector (F_in x 1)
            % @epsilon: Self-loop scaling factor

            % author: Anne Tumlin
            % date: 1/13/2026

            switch nargin

                case 6
                    name = varargin{1};
                    W_node = varargin{2};
                    b_node = varargin{3};
                    W_edge = varargin{4};
                    b_edge = varargin{5};
                    epsilon = varargin{6};

                    obj = obj.initWithWeights(name, W_node, b_node, W_edge, b_edge, epsilon);

                case 5
                    name = varargin{1};
                    W_node = varargin{2};
                    b_node = varargin{3};
                    W_edge = varargin{4};
                    b_edge = varargin{5};

                    obj = obj.initWithWeights(name, W_node, b_node, W_edge, b_edge, 0);

                case 4
                    W_node = varargin{1};
                    b_node = varargin{2};
                    W_edge = varargin{3};
                    b_edge = varargin{4};

                    obj = obj.initWithWeights('gine_layer', W_node, b_node, W_edge, b_edge, 0);

                case 0
                    obj.Name = 'gine_layer';
                    obj.InputSize = 0;
                    obj.OutputSize = 0;
                    obj.EdgeInputSize = 0;
                    obj.Weights = [];
                    obj.Bias = [];
                    obj.EdgeWeights = [];
                    obj.EdgeBias = [];
                    obj.Epsilon = 0;

                otherwise
                    error('Invalid number of inputs (should be 0, 4, 5, or 6)');
            end

        end

        % helper for constructor
        function obj = initWithWeights(obj, name, W_node, b_node, W_edge, b_edge, epsilon)
            if ~ischar(name) && ~isstring(name)
                error('Name must be a string or char array');
            end
            obj.Name = name;

            % Validate node weights
            if size(W_node, 2) ~= size(b_node, 1)
                error('Inconsistent node dimensions: size(W_node,2) must equal size(b_node,1)');
            end
            if size(b_node, 2) ~= 1
                error('Node bias vector should have one column');
            end

            % Validate edge weights - W_edge projects E_in to F_out
            if size(W_edge, 2) ~= size(W_node, 2)
                error('Edge output dimension must match node output dimension: size(W_edge,2) must equal size(W_node,2)');
            end
            if size(b_edge, 1) ~= size(W_edge, 2)
                error('Edge bias dimension must match edge output: size(b_edge,1) must equal size(W_edge,2)');
            end
            if size(b_edge, 2) ~= 1
                error('Edge bias vector should have one column');
            end

            obj.InputSize = size(W_node, 1);
            obj.OutputSize = size(W_node, 2);
            obj.EdgeInputSize = size(W_edge, 1);
            obj.Weights = W_node;
            obj.Bias = b_node;
            obj.EdgeWeights = W_edge;
            obj.EdgeBias = b_edge;
            obj.Epsilon = epsilon;
        end

        % evaluation method
        function Y = evaluate(obj, X, E, adj_list, edge_weights)
            % Forward pass through GINE layer (GNNV-compatible architecture)
            % @X: Input node features (N x F_in)
            % @E: Edge features (m x E_in)
            % @adj_list: Edge list (m x 2) with [src, dst] pairs
            % @edge_weights: Normalized adjacency weights per edge (m x 1)
            %                Optional - defaults to ones if not provided
            % @Y: Output node features (N x F_out)
            %
            % Architecture matches GNNV:
            %   1. If F_in != F_out: Transform nodes first with W
            %   2. Gather source features to edges
            %   3. Transform edge features to F_out and add
            %   4. ReLU
            %   5. Aggregate to destination nodes (weighted by edge_weights)
            %   6. Residual connection
            %   7. If F_in == F_out: Apply W at the end

            % author: Anne Tumlin
            % date: 1/13/2026

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
            edge_weights = edge_weights(:);  % Ensure column vector

            src_nodes = adj_list(:, 1);
            dst_nodes = adj_list(:, 2);

            % (1) Transform nodes first if dimensions change
            used_W = false;
            if obj.InputSize ~= obj.OutputSize
                H = X * obj.Weights;  % [N x F_out]
                used_W = true;
            else
                H = X;  % [N x F_out]
            end

            % (2) Gather source node features to edges
            H_src = H(src_nodes, :);  % [m x F_out]

            % (3) Transform edge features to F_out and add
            E_trans = E * obj.EdgeWeights;  % [m x F_out]
            if ~isempty(obj.EdgeBias)
                E_trans = E_trans + obj.EdgeBias';
            end
            edge_sum = H_src + E_trans;  % [m x F_out]

            % (4) ReLU
            edge_msg = max(0, edge_sum);  % [m x F_out]

            % (5) Aggregate messages to destination nodes (weighted)
            agg = zeros(numNodes, obj.OutputSize);
            for e = 1:numEdges
                agg(dst_nodes(e), :) = agg(dst_nodes(e), :) + edge_weights(e) * edge_msg(e, :);
            end

            % (6) Residual connection: H + aggregated
            combined = H + agg;  % [N x F_out]

            % (7) Apply W at end if not used earlier
            if ~used_W
                Y = combined * obj.Weights;  % [N x F_out]
            else
                Y = combined;
            end

            % Add bias
            if ~isempty(obj.Bias)
                Y = Y + obj.Bias';
            end

        end

        % main reachability analysis function
        function S = reach(varargin)
            % Reachability analysis for GINELayer
            % @in_set: Input GraphStar set
            % @E: Edge features (matrix for node-only, Star for edge perturbation)
            % @adj_list: Edge list (m x 2)
            % @method: 'approx-star', 'exact-star', 'approx-zono'
            % @option: 'single' or 'parallel'
            % @edge_weights: Edge weights for aggregation (m x 1) - optional

            % author: Anne Tumlin
            % date: 1/13/2026

            edge_weights = [];  % default

            switch nargin

                case 9
                    obj = varargin{1};
                    in_sets = varargin{2};
                    E = varargin{3};
                    adj_list = varargin{4};
                    method = varargin{5};
                    option = varargin{6};
                    % relaxFactor = varargin{7}; % not used
                    % lp_solver = varargin{8}; % not used
                    edge_weights = varargin{9};

                case 8
                    obj = varargin{1};
                    in_sets = varargin{2};
                    E = varargin{3};
                    adj_list = varargin{4};
                    method = varargin{5};
                    option = varargin{6};
                    % relaxFactor = varargin{7}; % not used
                    % lp_solver = varargin{8}; % not used

                case 7
                    obj = varargin{1};
                    in_sets = varargin{2};
                    E = varargin{3};
                    adj_list = varargin{4};
                    method = varargin{5};
                    option = varargin{6};
                    % relaxFactor = varargin{7}; % not used

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
                error('Zonotope reachability for GINELayer not yet implemented');
            else
                error('Unknown reachability method: %s', method);
            end

        end

    end


    methods % reachability methods

        % reachability for a single GraphStar input
        function gs_out = reach_star_single_input(obj, in_gs, E, adj_list, edge_weights)
            % Reachability through GINE layer for single GraphStar
            % @in_gs: Input GraphStar set
            % @E: Edge features (matrix for node-only mode)
            % @adj_list: Edge list (m x 2)
            % @edge_weights: Edge weights for aggregation (m x 1) - optional
            % @gs_out: Output GraphStar set

            % author: Anne Tumlin
            % date: 1/13/2026

            if nargin < 5
                edge_weights = [];
            end

            % Type detection for perturbation mode
            if isa(E, 'Star') || isa(E, 'ImageStar') || isa(E, 'GraphStar')
                % Edge perturbation mode - E is a Star/ImageStar/GraphStar
                gs_out = obj.reach_with_edge_perturbation(in_gs, E, adj_list, edge_weights);
            else
                % Node-only mode - E is constant matrix
                gs_out = obj.reach_node_only(in_gs, E, adj_list, edge_weights);
            end

        end

        % node-only perturbation reachability
        function gs_out = reach_node_only(obj, in_gs, E, adj_list, edge_weights)
            % Node-only reachability (edge features as constants)
            % GNNV-compatible architecture:
            %   1. If F_in != F_out: Transform nodes first with W
            %   2. Gather source features to edges
            %   3. Add edge encoding (constant)
            %   4. ReLU (using NNV's standard ReluLayer)
            %   5. Aggregate to nodes (weighted)
            %   6. Residual + final transform (if needed)
            %
            % @in_gs: Input GraphStar set
            % @E: Edge features (m x E_in) - constant matrix
            % @adj_list: Edge list (m x 2)
            % @edge_weights: Edge weights for aggregation (m x 1) - optional
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
            numEdges = size(adj_list, 1);

            % Default edge weights to 1 if not provided
            if nargin < 5 || isempty(edge_weights)
                edge_weights = ones(numEdges, 1);
            end
            edge_weights = edge_weights(:);

            src_nodes = adj_list(:, 1);
            dst_nodes = adj_list(:, 2);

            % (1) Transform nodes first if dimensions change
            used_W = false;
            if obj.InputSize ~= obj.OutputSize
                V_H = obj.map_features_star(in_gs.V, obj.Weights);  % [N x F_out x K]
                used_W = true;
            else
                V_H = in_gs.V;  % [N x F_out x K]
            end

            % (2) Gather source node features to edges
            V_edge = obj.gather_src_features_star(V_H, src_nodes);  % [m x F_out x K]

            % (3) Transform edge features to F_out (constant) and add to center
            E_trans = E * obj.EdgeWeights;  % [m x F_out]
            if ~isempty(obj.EdgeBias)
                E_trans = E_trans + obj.EdgeBias';
            end
            V_edge(:, :, 1) = V_edge(:, :, 1) + E_trans;

            % (4) Apply ReLU using NNV's standard ReluLayer (sound implementation)
            % Convert GraphStar V [m x F x K] to Star V [m*F x K]
            [m_edges, F_out, K_orig] = size(V_edge);
            V_flat = reshape(permute(V_edge, [2, 1, 3]), [m_edges * F_out, K_orig]);

            % Create Star for edge features
            edge_star = Star(V_flat, in_gs.C, in_gs.d, in_gs.pred_lb, in_gs.pred_ub);

            % Apply ReLU using standard NNV ReluLayer
            L = ReluLayer();
            edge_star_out = L.reach(edge_star, 'approx-star');

            % Handle cell array output (approx-star can return cell)
            if iscell(edge_star_out)
                edge_star_out = edge_star_out{1};
            end

            % Reshape back to [m x F x K']
            K_new = edge_star_out.nVar + 1;
            V_flat_out = edge_star_out.V;
            V_edge_relu = permute(reshape(V_flat_out, [F_out, m_edges, K_new]), [2, 1, 3]);

            % Get updated constraints
            C_new = edge_star_out.C;
            d_new = edge_star_out.d;
            pred_lb_new = edge_star_out.predicate_lb;
            pred_ub_new = edge_star_out.predicate_ub;

            % (5) Aggregate edge messages back to nodes (weighted)
            V_agg = obj.aggregate_edges_to_nodes_weighted_star(V_edge_relu, dst_nodes, numNodes, edge_weights);

            % (6) Residual connection: H + aggregated
            % Need to match K dimensions
            K_agg = size(V_agg, 3);
            K_H = size(V_H, 3);
            if K_agg > K_H
                V_H_expanded = zeros(numNodes, obj.OutputSize, K_agg, 'like', V_H);
                V_H_expanded(:, :, 1:K_H) = V_H;
                V_combined = V_H_expanded + V_agg;
            else
                V_combined = V_H + V_agg;
            end

            % (7) Apply W at end if not used earlier
            if ~used_W
                V_out = obj.map_features_star(V_combined, obj.Weights);
            else
                V_out = V_combined;
            end

            % Add bias to center only
            if ~isempty(obj.Bias)
                V_out(:, :, 1) = V_out(:, :, 1) + obj.Bias';
            end

            % Create output GraphStar
            gs_out = GraphStar(V_out, C_new, d_new, pred_lb_new, pred_ub_new);

        end

        % edge perturbation reachability
        function gs_out = reach_with_edge_perturbation(obj, in_gs, E_star, adj_list, edge_weights)
            % Edge perturbation reachability (both node and edge features perturbed)
            % @in_gs: Input GraphStar set for node features
            % @E_star: Edge features as Star, ImageStar, or GraphStar (m x E_in)
            % @adj_list: Edge list (m x 2)
            % @edge_weights: Edge weights for aggregation (m x 1) - optional
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
            numEdges = size(adj_list, 1);

            % Default edge weights to 1 if not provided
            if nargin < 5 || isempty(edge_weights)
                edge_weights = ones(numEdges, 1);
            end
            edge_weights = edge_weights(:);

            src_nodes = adj_list(:, 1);
            dst_nodes = adj_list(:, 2);

            % (1) Transform node features first if dimensions change (GNNV architecture)
            used_W = false;
            if obj.InputSize ~= obj.OutputSize
                V_H = obj.map_features_star(in_gs.V, obj.Weights);  % [N x F_out x K]
                used_W = true;
            else
                V_H = in_gs.V;  % [N x F_out x K]
            end

            % Extract edge Star properties
            if isa(E_star, 'ImageStar')
                % ImageStar has V with shape [m x E_in x 1 x K] or similar
                E_V = squeeze(E_star.V);  % [m x E_in x K] or [m x E_in]
                if ndims(E_V) == 2
                    E_V = reshape(E_V, [size(E_V, 1), size(E_V, 2), 1]);
                end
                E_C = E_star.C;
                E_d = E_star.d;
                E_pred_lb = E_star.pred_lb;
                E_pred_ub = E_star.pred_ub;
            elseif isa(E_star, 'GraphStar')
                % GraphStar has V with shape [m x E_in x K] - same as node GraphStar
                E_V = E_star.V;  % [m x E_in x K]
                E_C = E_star.C;
                E_d = E_star.d;
                E_pred_lb = E_star.pred_lb;
                E_pred_ub = E_star.pred_ub;
            else
                % Star - flatten and reshape
                E_V = reshape(E_star.V, [numEdges, obj.EdgeInputSize, E_star.nVar + 1]);
                E_C = E_star.C;
                E_d = E_star.d;
                E_pred_lb = E_star.predicate_lb;
                E_pred_ub = E_star.predicate_ub;
            end

            % (2) Gather source node features to edges (now in F_out space)
            % V_node_edge: [m x F_out x K_node]
            V_node_edge = obj.gather_src_features_star(V_H, src_nodes);

            % (3) Transform edge features through edge weights to F_out
            % E_trans_V: [m x F_out x K_edge]
            E_trans_V = obj.map_features_star(E_V, obj.EdgeWeights);
            if ~isempty(obj.EdgeBias)
                E_trans_V(:, :, 1) = E_trans_V(:, :, 1) + obj.EdgeBias';
            end

            % (4) Combine node and edge features via Minkowski sum
            % Need to pad both to same number of generators
            [V_node_padded, C_node_padded, d_node, pred_lb_node, pred_ub_node, ...
             V_edge_padded, C_edge_padded, d_edge, pred_lb_edge, pred_ub_edge] = ...
                obj.pad_to_same_K(V_node_edge, in_gs.C, in_gs.d, in_gs.pred_lb, in_gs.pred_ub, ...
                                  E_trans_V, E_C, E_d, E_pred_lb, E_pred_ub);

            % Minkowski sum: add centers, concatenate generators
            K_node = size(V_node_padded, 3);
            K_edge = size(V_edge_padded, 3);

            % Combined V: center is sum of centers, generators are concatenated
            % Now both are in F_out space
            V_combined = zeros(numEdges, obj.OutputSize, K_node + K_edge - 1, 'like', V_node_padded);
            V_combined(:, :, 1) = V_node_padded(:, :, 1) + V_edge_padded(:, :, 1);  % centers
            V_combined(:, :, 2:K_node) = V_node_padded(:, :, 2:end);  % node generators
            V_combined(:, :, K_node+1:end) = V_edge_padded(:, :, 2:end);  % edge generators

            % Combined constraints: block diagonal
            C_combined = blkdiag(C_node_padded, C_edge_padded);
            d_combined = [d_node; d_edge];
            pred_lb_combined = [pred_lb_node; pred_lb_edge];
            pred_ub_combined = [pred_ub_node; pred_ub_edge];

            % (5) Apply ReLU using NNV's standard ReluLayer (sound implementation)
            % Convert V [m x F x K] to Star V [m*F x K]
            [m_edges, F_out, K_orig] = size(V_combined);
            V_flat = reshape(permute(V_combined, [2, 1, 3]), [m_edges * F_out, K_orig]);

            % Create Star for combined edge features
            edge_star = Star(V_flat, C_combined, d_combined, pred_lb_combined, pred_ub_combined);

            % Apply ReLU using standard NNV ReluLayer
            L = ReluLayer();
            edge_star_out = L.reach(edge_star, 'approx-star');

            % Handle cell array output (approx-star can return cell)
            if iscell(edge_star_out)
                edge_star_out = edge_star_out{1};
            end

            % Reshape back to [m x F x K']
            K_new = edge_star_out.nVar + 1;
            V_flat_out = edge_star_out.V;
            V_edge_relu = permute(reshape(V_flat_out, [F_out, m_edges, K_new]), [2, 1, 3]);

            % Get updated constraints
            C_new = edge_star_out.C;
            d_new = edge_star_out.d;
            pred_lb_new = edge_star_out.predicate_lb;
            pred_ub_new = edge_star_out.predicate_ub;

            % (6) Aggregate edge messages back to nodes (weighted)
            V_agg = obj.aggregate_edges_to_nodes_weighted_star(V_edge_relu, dst_nodes, numNodes, edge_weights);

            % (7) Residual connection: V_H + aggregated
            % Need to expand V_H to match V_agg's generator count
            K_agg = size(V_agg, 3);
            K_H = size(V_H, 3);
            if K_agg > K_H
                V_H_expanded = zeros(numNodes, obj.OutputSize, K_agg, 'like', V_H);
                V_H_expanded(:, :, 1:K_H) = V_H;
                V_combined_final = V_H_expanded + V_agg;
            else
                V_combined_final = V_H + V_agg;
            end

            % (8) Apply W at end if not used earlier (GNNV architecture)
            if ~used_W
                V_out = obj.map_features_star(V_combined_final, obj.Weights);
            else
                V_out = V_combined_final;
            end

            % Add bias to center only
            if ~isempty(obj.Bias)
                V_out(:, :, 1) = V_out(:, :, 1) + obj.Bias';
            end

            % Create output GraphStar
            gs_out = GraphStar(V_out, C_new, d_new, pred_lb_new, pred_ub_new);

        end

        % reachability for multiple GraphStar inputs
        function S = reach_star_multipleInputs(obj, in_sets, E, adj_list, option, edge_weights)
            % Reachability for multiple GraphStar inputs
            % @in_sets: Array of GraphStar sets
            % @E: Edge features
            % @adj_list: Edge list
            % @option: 'single' or 'parallel'
            % @edge_weights: Edge weights for aggregation (m x 1) - optional

            % author: Anne Tumlin
            % date: 1/13/2026

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

        % Aggregate edge features back to nodes
        function V_node = aggregate_edges_to_nodes_star(~, V_edge, dst_nodes, numNodes)
            % V_edge: [m x F x K] -> V_node: [N x F x K]
            [m, F, K] = size(V_edge);
            V_node = zeros(numNodes, F, K, 'like', V_edge);
            for k = 1:K
                for e = 1:m
                    V_node(dst_nodes(e), :, k) = V_node(dst_nodes(e), :, k) + V_edge(e, :, k);
                end
            end
        end

        % Aggregate edge features back to nodes with edge weights
        function V_node = aggregate_edges_to_nodes_weighted_star(~, V_edge, dst_nodes, numNodes, edge_weights)
            % V_edge: [m x F x K] -> V_node: [N x F x K]
            % edge_weights: [m x 1] - weights for each edge
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
            % Pad both V matrices to have same K dimension
            % Returns padded versions with consistent dimensions

            K1 = size(V1, 3);
            K2 = size(V2, 3);

            if K1 == K2
                % Already same size
                V1_out = V1;
                V2_out = V2;
                C1_out = C1;
                C2_out = C2;
                d1_out = d1;
                d2_out = d2;
                lb1_out = lb1;
                ub1_out = ub1;
                lb2_out = lb2;
                ub2_out = ub2;
                return;
            end

            % Ensure constraints are properly sized
            if isempty(C1)
                C1 = zeros(0, K1 - 1);
            end
            if isempty(C2)
                C2 = zeros(0, K2 - 1);
            end
            if isempty(d1)
                d1 = zeros(0, 1);
            end
            if isempty(d2)
                d2 = zeros(0, 1);
            end
            if isempty(lb1)
                lb1 = -ones(K1 - 1, 1);
            end
            if isempty(ub1)
                ub1 = ones(K1 - 1, 1);
            end
            if isempty(lb2)
                lb2 = -ones(K2 - 1, 1);
            end
            if isempty(ub2)
                ub2 = ones(K2 - 1, 1);
            end

            % Actually pad to same K
            K_max = max(K1, K2);

            if K1 < K_max
                % Pad V1 with zeros
                V1_out = zeros(size(V1, 1), size(V1, 2), K_max, 'like', V1);
                V1_out(:, :, 1:K1) = V1;
                % Pad constraints with zeros for new predicates (they're unconstrained)
                C1_out = [C1, zeros(size(C1, 1), K_max - K1)];
                d1_out = d1;
                lb1_out = [lb1; zeros(K_max - K1, 1)];  % new preds bounded [0,0]
                ub1_out = [ub1; zeros(K_max - K1, 1)];
            else
                V1_out = V1;
                C1_out = C1;
                d1_out = d1;
                lb1_out = lb1;
                ub1_out = ub1;
            end

            if K2 < K_max
                % Pad V2 with zeros
                V2_out = zeros(size(V2, 1), size(V2, 2), K_max, 'like', V2);
                V2_out(:, :, 1:K2) = V2;
                % Pad constraints with zeros for new predicates
                C2_out = [C2, zeros(size(C2, 1), K_max - K2)];
                d2_out = d2;
                lb2_out = [lb2; zeros(K_max - K2, 1)];
                ub2_out = [ub2; zeros(K_max - K2, 1)];
            else
                V2_out = V2;
                C2_out = C2;
                d2_out = d2;
                lb2_out = lb2;
                ub2_out = ub2;
            end

        end

    end


    methods % helper methods

        function obj = toGPU(obj)
            % Move weights to GPU

            obj.Weights = gpuArray(obj.Weights);
            obj.EdgeWeights = gpuArray(obj.EdgeWeights);
            if ~isempty(obj.Bias)
                obj.Bias = gpuArray(obj.Bias);
            end
            if ~isempty(obj.EdgeBias)
                obj.EdgeBias = gpuArray(obj.EdgeBias);
            end
        end

        function obj = toCPU(obj)
            % Move weights to CPU

            obj.Weights = gather(obj.Weights);
            obj.EdgeWeights = gather(obj.EdgeWeights);
            if ~isempty(obj.Bias)
                obj.Bias = gather(obj.Bias);
            end
            if ~isempty(obj.EdgeBias)
                obj.EdgeBias = gather(obj.EdgeBias);
            end
        end

        function obj = changeParamsPrecision(obj, precision)
            % Change precision of parameters
            % @precision: 'single' or 'double'

            if strcmp(precision, 'double')
                obj.Weights = double(obj.Weights);
                obj.EdgeWeights = double(obj.EdgeWeights);
                if ~isempty(obj.Bias)
                    obj.Bias = double(obj.Bias);
                end
                if ~isempty(obj.EdgeBias)
                    obj.EdgeBias = double(obj.EdgeBias);
                end
            elseif strcmp(precision, 'single')
                obj.Weights = single(obj.Weights);
                obj.EdgeWeights = single(obj.EdgeWeights);
                if ~isempty(obj.Bias)
                    obj.Bias = single(obj.Bias);
                end
                if ~isempty(obj.EdgeBias)
                    obj.EdgeBias = single(obj.EdgeBias);
                end
            else
                error('Precision must be "single" or "double"');
            end
        end

    end

end
