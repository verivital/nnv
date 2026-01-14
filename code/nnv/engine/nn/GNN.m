classdef GNN < handle
    % GNN class encodes Graph Neural Networks for verification in NNV
    %   This class wraps GNN layers (GCNLayer, GINELayer) and manages
    %   graph structure separately from layer weights, enabling weight
    %   reuse across different graphs.
    %
    %   Unlike standard NN layers, GNN layers require graph structure
    %   (adjacency, edge list, edge features) which this wrapper provides.
    %
    % Main references:
    % 1) Graph Neural Networks: A Review of Methods and Applications
    %    https://arxiv.org/abs/1812.08434
    % 2) NNV verification tool documentation
    %    https://github.com/verivital/nnv
    %
    % Author: Anne Tumlin
    % Date: 01/13/2026

    properties
        Name = 'gnn';          % Name of the network
        Layers = {};           % Cell array of GNN layers
        numLayers = 0;         % Number of layers
        InputSize = 0;         % Input feature dimension
        OutputSize = 0;        % Output feature dimension

        % GNN-specific: Graph structure
        A_norm = [];           % Normalized adjacency matrix (for GCN layers)
        adj_list = [];         % Edge list [src, dst] (for GINE layers)
        E = [];                % Edge features (for GINE layers)
        edge_weights = [];     % Edge weights for aggregation (for GINE layers)

        % Reachability options (mirrors NN.m)
        reachMethod = 'approx-star';
        relaxFactor = 0;
        reachOption = [];
        numCores = 1;
        reachSet = {};         % Reachable set for each layer
        reachTime = [];        % Computation time for each layer
        dis_opt = [];          % Display option
        lp_solver = 'linprog'; % LP solver for reachability
    end


    methods % Main methods (constructor, evaluation, reach)

        function obj = GNN(varargin)
            % GNN Constructor
            %
            % Syntax:
            %   gnn = GNN()                                        % Empty GNN
            %   gnn = GNN(layers)                                  % Layers only
            %   gnn = GNN(layers, A_norm)                          % GCN-only network
            %   gnn = GNN(layers, A_norm, adj_list, E)             % GINE without weights
            %   gnn = GNN(layers, A_norm, adj_list, E, edge_wts)   % GINE with weights
            %   gnn = GNN(layers, A_norm, adj_list, E, edge_wts, name)
            %
            % Inputs:
            %   layers     - Cell array of GNN layers
            %   A_norm     - Normalized adjacency matrix (N x N)
            %   adj_list   - Edge list (m x 2) with [src, dst] pairs
            %   E          - Edge features (m x E_in)
            %   edge_wts   - Edge weights for aggregation (m x 1)
            %   name       - Network name (string)

            switch nargin
                case 0
                    % Empty constructor
                    return;

                case 1
                    % Layers only
                    obj.Layers = varargin{1};

                case 2
                    % Layers + A_norm (GCN-only)
                    obj.Layers = varargin{1};
                    obj.A_norm = varargin{2};

                case 4
                    % Layers + A_norm + adj_list + E (no edge weights)
                    obj.Layers = varargin{1};
                    obj.A_norm = varargin{2};
                    obj.adj_list = varargin{3};
                    obj.E = varargin{4};

                case 5
                    % Layers + A_norm + adj_list + E + edge_weights
                    obj.Layers = varargin{1};
                    obj.A_norm = varargin{2};
                    obj.adj_list = varargin{3};
                    obj.E = varargin{4};
                    obj.edge_weights = varargin{5};

                case 6
                    % Layers + A_norm + adj_list + E + edge_weights + name
                    obj.Layers = varargin{1};
                    obj.A_norm = varargin{2};
                    obj.adj_list = varargin{3};
                    obj.E = varargin{4};
                    obj.edge_weights = varargin{5};
                    obj.Name = varargin{6};

                otherwise
                    error('Invalid number of inputs. Expected 0, 1, 2, 4, 5, or 6 arguments.');
            end

            % Update layer count
            obj.numLayers = length(obj.Layers);

            % Determine input/output sizes from layers if available
            if obj.numLayers > 0
                firstLayer = obj.Layers{1};
                if isprop(firstLayer, 'InputSize')
                    obj.InputSize = firstLayer.InputSize;
                end
                lastLayer = obj.Layers{end};
                if isprop(lastLayer, 'OutputSize')
                    obj.OutputSize = lastLayer.OutputSize;
                end
            end
        end


        function Y = evaluate(obj, X, E_sample)
            % Evaluate GNN given input node features
            %
            % Syntax:
            %   Y = gnn.evaluate(X)
            %   Y = gnn.evaluate(X, E_sample)  % For edge perturbation sampling
            %
            % Inputs:
            %   X        - Node features (N x F_in)
            %   E_sample - (Optional) Edge features (m x E_in) for sampling
            %              If not provided, uses stored obj.E
            %
            % Outputs:
            %   Y - Output features (N x F_out)

            % Use provided edge features or stored ones
            if nargin < 3 || isempty(E_sample)
                E_eval = obj.E;
                % If E is a Star/ImageStar/GraphStar (edge perturbation mode),
                % extract the center for evaluation
                if isa(E_eval, 'GraphStar')
                    E_eval = E_eval.V(:, :, 1);  % Center of GraphStar
                elseif isa(E_eval, 'Star') || isa(E_eval, 'ImageStar')
                    E_eval = E_eval.V(:, 1);  % Center of Star
                    E_eval = reshape(E_eval, size(obj.adj_list, 1), []);
                end
            else
                E_eval = E_sample;
            end

            Y = X;
            for i = 1:obj.numLayers
                layer = obj.Layers{i};

                if isa(layer, 'GCNLayer')
                    Y = layer.evaluate(Y, obj.A_norm);
                elseif isa(layer, 'GINELayer')
                    Y = layer.evaluate(Y, E_eval, obj.adj_list, obj.edge_weights);
                else
                    % Standard layer (ReLU, etc.)
                    Y = layer.evaluate(Y);
                end
            end
        end


        function outputSet = reach(obj, inputSet, reachOptions)
            % Compute reachable set through GNN
            %
            % Syntax:
            %   outputSet = gnn.reach(inputSet)
            %   outputSet = gnn.reach(inputSet, reachOptions)
            %
            % Inputs:
            %   inputSet     - GraphStar input set
            %   reachOptions - struct with fields:
            %       reachMethod  - 'approx-star' (default)
            %       numCores     - Number of cores (default: 1)
            %       relaxFactor  - Relaxation factor (default: 0)
            %       dis_opt      - Display option (default: [])
            %       lp_solver    - 'linprog' or 'glpk' (default: 'linprog')
            %
            % Outputs:
            %   outputSet - GraphStar output set

            % Validate input set type
            if ~isa(inputSet, 'GraphStar')
                error('Input must be a GraphStar. Got: %s', class(inputSet));
            end

            % Parse reachability options
            if exist('reachOptions', 'var') && isstruct(reachOptions)
                if isfield(reachOptions, 'reachMethod')
                    obj.reachMethod = reachOptions.reachMethod;
                end
                if isfield(reachOptions, 'numCores')
                    obj.numCores = reachOptions.numCores;
                end
                if isfield(reachOptions, 'relaxFactor')
                    obj.relaxFactor = reachOptions.relaxFactor;
                end
                if isfield(reachOptions, 'dis_opt')
                    obj.dis_opt = reachOptions.dis_opt;
                end
                if isfield(reachOptions, 'lp_solver')
                    obj.lp_solver = reachOptions.lp_solver;
                end
            end

            % Start parallel pool if needed
            if obj.numCores > 1
                obj.start_pool();
                obj.reachOption = 'parallel';
            else
                obj.reachOption = 'single';
            end

            % Display option
            if strcmp(obj.dis_opt, 'display')
                fprintf('\nPerform reachability analysis for GNN: %s\n', obj.Name);
            end

            % Initialize storage
            obj.reachSet = cell(1, obj.numLayers);
            obj.reachTime = zeros(1, obj.numLayers);

            % Layer-by-layer reachability
            rs = inputSet;
            for i = 1:obj.numLayers
                layer = obj.Layers{i};

                if strcmp(obj.dis_opt, 'display')
                    fprintf('  Layer %d (%s)...', i, layer.Name);
                end

                t = tic;

                if isa(layer, 'GCNLayer')
                    rs = layer.reach(rs, obj.A_norm, obj.reachMethod, obj.reachOption);
                elseif isa(layer, 'GINELayer')
                    rs = layer.reach(rs, obj.E, obj.adj_list, obj.reachMethod, obj.reachOption, obj.relaxFactor, obj.lp_solver, obj.edge_weights);
                elseif isa(layer, 'ReluLayer') || isa(layer, 'ActivationFunctionLayer')
                    % Activation layers need Star, not GraphStar
                    % Convert GraphStar to Star, apply activation, convert back
                    numNodes = rs.numNodes;
                    numFeatures = rs.numFeatures;
                    S = rs.toStar();
                    S_out = layer.reach(S, obj.reachMethod, obj.reachOption, ...
                                       obj.relaxFactor, obj.dis_opt, obj.lp_solver);
                    % Handle cell array output from ReLU (can produce multiple sets)
                    if iscell(S_out)
                        S_out = S_out{1};  % Take first set for approx methods
                    end
                    % Convert back to GraphStar
                    new_V = reshape(S_out.V, [numNodes, numFeatures, S_out.nVar + 1]);
                    if ~isempty(S_out.state_lb) && ~isempty(S_out.state_ub)
                        nf_lb = reshape(S_out.state_lb, [numNodes, numFeatures]);
                        nf_ub = reshape(S_out.state_ub, [numNodes, numFeatures]);
                        rs = GraphStar(new_V, S_out.C, S_out.d, S_out.predicate_lb, S_out.predicate_ub, nf_lb, nf_ub);
                    else
                        rs = GraphStar(new_V, S_out.C, S_out.d, S_out.predicate_lb, S_out.predicate_ub);
                    end
                else
                    % Standard layer
                    rs = layer.reach(rs, obj.reachMethod, obj.reachOption, ...
                                    obj.relaxFactor, obj.dis_opt, obj.lp_solver);
                end

                obj.reachTime(i) = toc(t);
                obj.reachSet{i} = rs;

                if strcmp(obj.dis_opt, 'display')
                    fprintf(' done (%.4f sec)\n', obj.reachTime(i));
                end
            end

            outputSet = rs;
        end


        function obj = setGraph(obj, A_norm, adj_list, E)
            % Update graph structure (enables weight reuse on different graphs)
            %
            % Syntax:
            %   gnn.setGraph(A_norm)
            %   gnn.setGraph(A_norm, adj_list)
            %   gnn.setGraph(A_norm, adj_list, E)
            %
            % Inputs:
            %   A_norm   - Normalized adjacency matrix
            %   adj_list - Edge list [src, dst]
            %   E        - Edge features

            obj.A_norm = A_norm;
            if nargin >= 3
                obj.adj_list = adj_list;
            end
            if nargin >= 4
                obj.E = E;
            end
        end

    end


    methods % Helper methods

        function obj = params2gpu(obj)
            % Move all layer weights to GPU
            for i = 1:obj.numLayers
                if ismethod(obj.Layers{i}, 'toGPU')
                    obj.Layers{i} = obj.Layers{i}.toGPU();
                end
            end
        end


        function obj = changeParamsPrecision(obj, precision)
            % Change precision of all layer parameters
            %
            % Inputs:
            %   precision - 'single' or 'double'

            for i = 1:obj.numLayers
                if ismethod(obj.Layers{i}, 'changeParamsPrecision')
                    obj.Layers{i} = obj.Layers{i}.changeParamsPrecision(precision);
                end
            end
        end


        function start_pool(obj)
            % Start parallel pool for multi-core computation

            if obj.numCores > 1
                poolobj = gcp('nocreate');
                if isempty(poolobj)
                    parpool('local', obj.numCores);
                else
                    if poolobj.NumWorkers ~= obj.numCores
                        delete(poolobj);
                        parpool('local', obj.numCores);
                    end
                end
            end
        end


        function info = getInfo(obj)
            % Get summary information about the GNN
            %
            % Outputs:
            %   info - struct with network information

            info = struct();
            info.Name = obj.Name;
            info.numLayers = obj.numLayers;
            info.InputSize = obj.InputSize;
            info.OutputSize = obj.OutputSize;
            info.hasAdjacency = ~isempty(obj.A_norm);
            info.hasEdgeList = ~isempty(obj.adj_list);
            info.hasEdgeFeatures = ~isempty(obj.E);

            if ~isempty(obj.A_norm)
                info.numNodes = size(obj.A_norm, 1);
            end
            if ~isempty(obj.adj_list)
                info.numEdges = size(obj.adj_list, 1);
            end

            % Layer types
            info.layerTypes = cell(1, obj.numLayers);
            for i = 1:obj.numLayers
                info.layerTypes{i} = class(obj.Layers{i});
            end
        end


        function disp(obj)
            % Display GNN information

            fprintf('GNN: %s\n', obj.Name);
            fprintf('  Layers: %d\n', obj.numLayers);
            fprintf('  Input size: %d\n', obj.InputSize);
            fprintf('  Output size: %d\n', obj.OutputSize);

            if ~isempty(obj.A_norm)
                fprintf('  Nodes: %d\n', size(obj.A_norm, 1));
            end
            if ~isempty(obj.adj_list)
                fprintf('  Edges: %d\n', size(obj.adj_list, 1));
            end

            fprintf('  Layer types:\n');
            for i = 1:obj.numLayers
                fprintf('    %d: %s (%s)\n', i, obj.Layers{i}.Name, class(obj.Layers{i}));
            end
        end

    end

end
