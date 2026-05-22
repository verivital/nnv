function [gnn, test_data, norm_stats] = gnn2nnv(mat_path)
% gnn2nnv - Import a Python-trained GNN model into NNV
%
% Automatically detects model type (GCN, GINE-linear, GINE-full) from
% the .mat file structure and constructs the correct GNN object.
%
% Syntax:
%   gnn = gnn2nnv(mat_path)
%   [gnn, test_data] = gnn2nnv(mat_path)
%   [gnn, test_data, norm_stats] = gnn2nnv(mat_path)
%
% Inputs:
%   mat_path - Path to .mat file exported by gnn_training Python project
%
% Outputs:
%   gnn        - GNN object ready for evaluate() and reach()
%   test_data  - struct with fields: X, Y, E (if applicable), adj_list
%   norm_stats - struct with normalization parameters (if stored)
%
% Supported model types (auto-detected from .mat fields):
%   'gcn'        - GCN layers (uses ANorm_g adjacency)
%   'sage'       - SAGEConv layers (uses A_adj binary adjacency)
%   'gine_linear'- Simplified GINE with linear projections (GINELayer)
%   'gine_conv'  - Full GINE with MLPs (GINEConvLayer)
%   'gine_pretrain' - GIN+E from Hu et al. ICLR 2020 (HuGINEConvLayer)
%
% Author: Anne Tumlin
% Date: 02/11/2026

    %% Validate input
    if nargin < 1
        error('gnn2nnv:noInput', 'Path to .mat file is required.');
    end
    if ~isfile(mat_path)
        error('gnn2nnv:fileNotFound', 'File not found: %s', mat_path);
    end

    %% Load .mat file
    model = load(mat_path);

    %% Detect model type
    model_type = detect_model_type(model);
    fprintf('gnn2nnv: Detected model type: %s\n', model_type);

    %% Build GNN based on detected type
    switch model_type
        case 'gcn'
            gnn = build_gcn(model);
        case 'sage'
            gnn = build_sage(model);
        case 'gine_linear'
            gnn = build_gine_linear(model);
        case 'gine_conv'
            gnn = build_gine_conv(model);
        case 'gine_pretrain'
            gnn = build_gine_pretrain(model);
        otherwise
            error('gnn2nnv:unknownType', 'Unknown model type: %s', model_type);
    end

    %% Cross-validate against Python predictions if available
    if isfield(model, 'python_predictions')
        validate_predictions(gnn, model, model_type);
    end

    %% Extract test data
    test_data = extract_test_data(model, model_type);

    %% Extract normalization stats
    norm_stats = extract_norm_stats(model);

    fprintf('gnn2nnv: Import complete. GNN has %d layers.\n', gnn.numLayers);
end


%% ======== Model Type Detection ========

function model_type = detect_model_type(model)
% Detect GNN model type from .mat file fields

    % Explicit model_type field takes priority
    if isfield(model, 'model_type')
        model_type = char(model.model_type);
        return;
    end

    params = [];
    if isfield(model, 'best_params')
        params = model.best_params;
    end

    if isempty(params)
        error('gnn2nnv:noParams', 'No best_params found in .mat file.');
    end

    % Check for Hu et al. GIN+E (has conv.edge_proj, not edge_linear)
    if has_conv_edge_proj_fields(params)
        model_type = 'gine_pretrain';
        return;
    end

    % Check for full GINE (has conv.mlp fields with edge_linear)
    if has_conv_mlp_fields(params)
        model_type = 'gine_conv';
        return;
    end

    % Check for simplified GINE (has edge fields but no MLP)
    if has_edge_fields(params)
        model_type = 'gine_linear';
        return;
    end

    % Check for SAGEConv (has sage{i} fields with NodeWeights)
    if has_sage_fields(params)
        model_type = 'sage';
        return;
    end

    % Check for GCN (has mult fields and ANorm_g)
    if has_mult_fields(params) && isfield(model, 'ANorm_g')
        model_type = 'gcn';
        return;
    end

    error('gnn2nnv:unrecognized', ...
        'Could not detect model type. Ensure .mat has model_type field or standard weight fields.');
end


function tf = has_mult_fields(params)
    tf = isfield(params, 'mult1') && isfield(params.mult1, 'Weights');
end


function tf = has_edge_fields(params)
    tf = isfield(params, 'edge1') && isfield(params.edge1, 'Weights');
end


function tf = has_conv_edge_proj_fields(params)
    tf = isfield(params, 'conv1') && isfield(params.conv1, 'edge_proj');
end


function tf = has_conv_mlp_fields(params)
    tf = isfield(params, 'conv1') && isfield(params.conv1, 'mlp1');
end


function tf = has_sage_fields(params)
    tf = isfield(params, 'sage1') && isfield(params.sage1, 'NodeWeights');
end


%% ======== GCN Builder ========

function gnn = build_gcn(model)
% Build GNN with GCN layers from .mat file
%
% Expected fields:
%   model.best_params.mult{i}.Weights  (F_in x F_out)
%   model.best_params.mult{i}.Bias     (F_out x 1) [optional]
%   model.ANorm_g                      (N x N)

    params = model.best_params;
    A_norm = double(model.ANorm_g);

    % Count layers
    num_layers = count_indexed_fields(params, 'mult');
    fprintf('gnn2nnv: Building GCN with %d conv layers\n', num_layers);

    % Determine if ReLU activations are interleaved
    has_relu = true;
    if isfield(model, 'activations')
        has_relu = ~strcmp(model.activations, 'none');
    end

    layers = {};
    for i = 1:num_layers
        field = sprintf('mult%d', i);
        W = double(gather_if_gpu(params.(field).Weights));
        if isfield(params.(field), 'Bias') && ~isempty(params.(field).Bias)
            b = double(gather_if_gpu(params.(field).Bias));
        else
            b = zeros(size(W, 2), 1);
        end

        name = sprintf('gcn%d', i);
        layers{end+1} = GCNLayer(name, W, b); %#ok<AGROW>
        fprintf('  Layer %d: GCNLayer %d -> %d\n', i, size(W, 1), size(W, 2));

        % Add ReLU after each conv layer (standard GCN architecture)
        if has_relu
            layers{end+1} = ReluLayer(); %#ok<AGROW>
        end
    end

    gnn = GNN(layers, A_norm);
end


%% ======== SAGEConv Builder ========

function gnn = build_sage(model)
% Build GNN with SAGEConv layers from .mat file
%
% Expected fields:
%   model.best_params.sage{i}.NodeWeights  (F_in x F_out)
%   model.best_params.sage{i}.EdgeWeights  (F_in x F_out)
%   model.best_params.sage{i}.Bias         (F_out x 1)
%   model.A_adj                            (N x N binary adjacency)

    params = model.best_params;
    A_adj = double(model.A_adj);

    % Count layers
    num_layers = count_indexed_fields(params, 'sage');
    fprintf('gnn2nnv: Building SAGEConv with %d conv layers\n', num_layers);

    % Determine if ReLU activations are interleaved
    has_relu = true;
    if isfield(model, 'activations')
        has_relu = ~strcmp(model.activations, 'none');
    end

    layers = {};
    for i = 1:num_layers
        field = sprintf('sage%d', i);
        W_node = double(gather_if_gpu(params.(field).NodeWeights));
        W_edge = double(gather_if_gpu(params.(field).EdgeWeights));
        if isfield(params.(field), 'Bias') && ~isempty(params.(field).Bias)
            b = double(gather_if_gpu(params.(field).Bias));
        else
            b = zeros(size(W_node, 2), 1);
        end

        name = sprintf('sage%d', i);
        layers{end+1} = SAGEConvLayer(name, W_node, W_edge, b); %#ok<AGROW>
        fprintf('  Layer %d: SAGEConvLayer %d -> %d\n', i, size(W_node, 1), size(W_node, 2));

        % Add ReLU after each conv layer
        if has_relu
            layers{end+1} = ReluLayer(); %#ok<AGROW>
        end
    end

    % A_adj stored in GNN's A_norm field (GNN.m passes it to SAGEConvLayer)
    gnn = GNN(layers, A_adj);
end


%% ======== Simplified GINE Builder ========

function gnn = build_gine_linear(model)
% Build GNN with GINELayer (simplified linear variant)
%
% Expected fields:
%   model.best_params.mult{i}.Weights  (F_in x F_out) node weights
%   model.best_params.edge{i}.Weights  (E_in x F_out) edge weights
%   model.src, model.dst               edge list (1-indexed)
%   model.E_edge                       edge features
%   model.a                            edge weights for aggregation

    params = model.best_params;

    % Extract graph structure
    src = double(model.src);
    dst = double(model.dst);
    adj_list = [src, dst];
    E = double(model.E_edge);
    edge_weights = [];
    if isfield(model, 'a')
        edge_weights = double(model.a);
    end

    % Count layers
    num_layers = count_indexed_fields(params, 'mult');
    fprintf('gnn2nnv: Building simplified GINE with %d layers\n', num_layers);

    layers = {};
    for i = 1:num_layers
        node_field = sprintf('mult%d', i);
        edge_field = sprintf('edge%d', i);

        W_node = double(gather_if_gpu(params.(node_field).Weights));
        if isfield(params.(node_field), 'Bias') && ~isempty(params.(node_field).Bias)
            b_node = double(gather_if_gpu(params.(node_field).Bias));
        else
            b_node = zeros(size(W_node, 2), 1);
        end

        W_edge = double(gather_if_gpu(params.(edge_field).Weights));
        if isfield(params.(edge_field), 'Bias') && ~isempty(params.(edge_field).Bias)
            b_edge = double(gather_if_gpu(params.(edge_field).Bias));
        else
            b_edge = zeros(size(W_edge, 2), 1);
        end

        name = sprintf('gine%d', i);
        layers{end+1} = GINELayer(name, W_node, b_node, W_edge, b_edge); %#ok<AGROW>
        fprintf('  Layer %d: GINELayer node %d->%d, edge %d->%d\n', ...
            i, size(W_node, 1), size(W_node, 2), size(W_edge, 1), size(W_edge, 2));
    end

    gnn = GNN(layers, [], adj_list, E, edge_weights);
end


%% ======== Full GINE (GINEConv) Builder ========

function gnn = build_gine_conv(model)
% Build GNN with GINEConvLayer (full GINE with MLPs)
%
% Expected fields:
%   model.best_params.conv{i}.mlp1.Weights     (F_in x hidden)
%   model.best_params.conv{i}.mlp1.Bias        (hidden x 1)
%   model.best_params.conv{i}.mlp2.Weights     (hidden x F_out)
%   model.best_params.conv{i}.mlp2.Bias        (F_out x 1)
%   model.best_params.conv{i}.edge_linear.Weights  (E_in x F_in)
%   model.best_params.conv{i}.edge_linear.Bias     (F_in x 1)
%   model.best_params.conv{i}.eps              (scalar, optional)
%   model.src, model.dst, model.E_edge, model.a

    params = model.best_params;

    % Extract graph structure
    src = double(model.src);
    dst = double(model.dst);
    adj_list = [src, dst];
    E = double(model.E_edge);
    edge_weights = [];
    if isfield(model, 'a')
        edge_weights = double(model.a);
    end

    % Count layers
    num_layers = count_indexed_fields(params, 'conv');
    fprintf('gnn2nnv: Building full GINE with %d conv layers\n', num_layers);

    % Determine if inter-layer ReLU activations are used
    has_relu = true;
    if isfield(model, 'activations')
        has_relu = ~strcmp(model.activations, 'none');
    end

    layers = {};
    for i = 1:num_layers
        conv_field = sprintf('conv%d', i);
        conv = params.(conv_field);

        % MLP layer 1
        W1 = double(gather_if_gpu(conv.mlp1.Weights));
        b1 = get_bias(conv.mlp1, size(W1, 2));

        % MLP layer 2
        W2 = double(gather_if_gpu(conv.mlp2.Weights));
        b2 = get_bias(conv.mlp2, size(W2, 2));

        % Edge projection
        W_edge = double(gather_if_gpu(conv.edge_linear.Weights));
        b_edge = get_bias(conv.edge_linear, size(W_edge, 2));

        % Epsilon (optional)
        eps_val = 0;
        if isfield(conv, 'eps')
            eps_val = double(conv.eps);
        end

        name = sprintf('gine_conv%d', i);
        layers{end+1} = GINEConvLayer(name, W1, b1, W2, b2, W_edge, b_edge, eps_val); %#ok<AGROW>
        fprintf('  Layer %d: GINEConvLayer F_in=%d, hidden=%d, F_out=%d, E_in=%d, eps=%.4f\n', ...
            i, size(W1, 1), size(W1, 2), size(W2, 2), size(W_edge, 1), eps_val);

        % Add inter-layer ReLU (matches Python: F.relu between conv layers)
        if has_relu && i < num_layers
            layers{end+1} = ReluLayer(); %#ok<AGROW>
        end
    end

    gnn = GNN(layers, [], adj_list, E, edge_weights);
end


%% ======== Hu et al. GIN+E (HuGINEConvLayer) Builder ========

function gnn = build_gine_pretrain(model)
% Build GNN with HuGINEConvLayer (Hu et al. ICLR 2020)
%
% Expected fields:
%   model.best_params.conv{i}.mlp1.Weights       (F_in x hidden)
%   model.best_params.conv{i}.mlp1.Bias           (hidden x 1)
%   model.best_params.conv{i}.mlp2.Weights        (hidden x F_out)
%   model.best_params.conv{i}.mlp2.Bias           (F_out x 1)
%   model.best_params.conv{i}.edge_proj.Weights   (E_in x F_in)
%   model.best_params.conv{i}.edge_proj.Bias      (F_in x 1)
%   model.src, model.dst, model.E_edge
%
% Self-loops should already be included in the edge list (added by
% the Python exporter). Edge features for self-loops are zeros.

    params = model.best_params;

    % Extract graph structure
    src = double(model.src);
    dst = double(model.dst);
    adj_list = [src, dst];
    E = double(model.E_edge);
    edge_weights = [];
    if isfield(model, 'a')
        edge_weights = double(model.a);
    end

    % Count layers
    num_layers = count_indexed_fields(params, 'conv');
    fprintf('gnn2nnv: Building Hu et al. GIN+E with %d conv layers\n', num_layers);

    layers = {};
    for i = 1:num_layers
        conv_field = sprintf('conv%d', i);
        conv = params.(conv_field);

        % MLP layer 1
        W1 = double(gather_if_gpu(conv.mlp1.Weights));
        b1 = get_bias(conv.mlp1, size(W1, 2));

        % MLP layer 2
        W2 = double(gather_if_gpu(conv.mlp2.Weights));
        b2 = get_bias(conv.mlp2, size(W2, 2));

        % Edge projection
        W_edge = double(gather_if_gpu(conv.edge_proj.Weights));
        b_edge = get_bias(conv.edge_proj, size(W_edge, 2));

        name = sprintf('hu_gine%d', i);
        layers{end+1} = HuGINEConvLayer(name, W1, b1, W2, b2, W_edge, b_edge); %#ok<AGROW>
        fprintf('  Layer %d: HuGINEConvLayer F_in=%d, hidden=%d, F_out=%d, E_in=%d\n', ...
            i, size(W1, 1), size(W1, 2), size(W2, 2), size(W_edge, 1));

        % Add inter-layer ReLU (outer ReLU for intermediate layers)
        if i < num_layers
            layers{end+1} = ReluLayer(); %#ok<AGROW>
        end
    end

    gnn = GNN(layers, [], adj_list, E, edge_weights);
end


%% ======== Validation ========

function validate_predictions(gnn, model, model_type)
% Cross-validate NNV evaluate() against stored Python predictions

    fprintf('gnn2nnv: Validating against Python predictions...\n');

    try
        % Get test input (graph 1 — validation is single-graph)
        if isfield(model, 'X_test_g')
            if iscell(model.X_test_g)
                X = double(model.X_test_g{1});
            else
                X = double(model.X_test_g);
            end
        else
            fprintf('  No test input found, skipping validation.\n');
            return;
        end

        if ~isfield(model, 'python_predictions')
            fprintf('  No python_predictions field, skipping validation.\n');
            return;
        end

        % Predictions may be stored per-graph as a cell array (mirroring
        % X_test_g) or as a single numeric matrix. Match the input.
        if iscell(model.python_predictions)
            py_pred = double(model.python_predictions{1});
        else
            py_pred = double(model.python_predictions);
        end
        nnv_pred = gnn.evaluate(X);

        max_diff = max(abs(nnv_pred(:) - py_pred(:)));
        fprintf('  Max difference (NNV vs Python): %.2e\n', max_diff);

        if max_diff > 1e-5
            warning('gnn2nnv:predMismatch', ...
                'NNV predictions differ from Python by %.2e (threshold: 1e-5).', max_diff);
        else
            fprintf('  Validation PASSED.\n');
        end
    catch ME
        warning('gnn2nnv:validationFailed', ...
            'Could not validate predictions: %s', ME.message);
    end
end


%% ======== Test Data Extraction ========

function test_data = extract_test_data(model, model_type)
% Extract test data from .mat file
% Supports multi-graph exports: X_all{i}, Y_all{i} for all graphs
% test_data.X and test_data.Y are graph 1 (backward compatible)

    test_data = struct();

    % Node features
    if isfield(model, 'X_test_g')
        if iscell(model.X_test_g)
            test_data.num_graphs = numel(model.X_test_g);
            test_data.X = double(model.X_test_g{1});
            % Store all graphs for multi-graph verification
            if test_data.num_graphs > 1
                test_data.X_all = cell(test_data.num_graphs, 1);
                for i = 1:test_data.num_graphs
                    test_data.X_all{i} = double(model.X_test_g{i});
                end
            else
                test_data.X_all = {test_data.X};
            end
        else
            test_data.num_graphs = 1;
            test_data.X = double(model.X_test_g);
            test_data.X_all = {test_data.X};
        end
    end

    % Labels
    if isfield(model, 'Y_test_g')
        if iscell(model.Y_test_g)
            test_data.Y = double(model.Y_test_g{1});
            if numel(model.Y_test_g) > 1
                test_data.Y_all = cell(numel(model.Y_test_g), 1);
                for i = 1:numel(model.Y_test_g)
                    test_data.Y_all{i} = double(model.Y_test_g{i});
                end
            else
                test_data.Y_all = {test_data.Y};
            end
        else
            test_data.Y = double(model.Y_test_g);
            test_data.Y_all = {test_data.Y};
        end
    end

    % Edge features and graph structure (for GINE models)
    if strcmp(model_type, 'gine_linear') || strcmp(model_type, 'gine_conv') || strcmp(model_type, 'gine_pretrain')
        if isfield(model, 'E_edge')
            test_data.E = double(model.E_edge);
        end
        if isfield(model, 'src') && isfield(model, 'dst')
            test_data.adj_list = [double(model.src), double(model.dst)];
        end
    end

    % Adjacency (for GCN)
    if strcmp(model_type, 'gcn') && isfield(model, 'ANorm_g')
        test_data.A_norm = double(model.ANorm_g);
    end

    % Adjacency (for SAGEConv)
    if strcmp(model_type, 'sage') && isfield(model, 'A_adj')
        test_data.A_adj = double(model.A_adj);
    end
end


%% ======== Normalization Stats Extraction ========

function norm_stats = extract_norm_stats(model)
% Extract normalization parameters if stored

    norm_stats = struct();

    norm_fields = {'X_mean', 'X_std', 'Y_mean', 'Y_std', ...
                   'X_max', 'Y_max', 'X_min', 'Y_min', ...
                   'x_mean', 'x_std', 'y_mean', 'y_std', ...
                   'norm_mean', 'norm_std', 'scaler_mean', 'scaler_std'};

    for i = 1:length(norm_fields)
        f = norm_fields{i};
        if isfield(model, f)
            norm_stats.(f) = double(model.(f));
        end
    end
end


%% ======== Helper Functions ========

function n = count_indexed_fields(params, prefix)
% Count sequential numbered fields (e.g., mult1, mult2, mult3 -> 3)
    n = 0;
    while isfield(params, sprintf('%s%d', prefix, n + 1))
        n = n + 1;
    end
end


function val = gather_if_gpu(val)
% Gather GPU array to CPU if needed
    if isa(val, 'gpuArray')
        val = gather(val);
    end
end


function b = get_bias(layer_struct, expected_size)
% Extract bias from struct, defaulting to zeros
    if isfield(layer_struct, 'Bias') && ~isempty(layer_struct.Bias)
        b = double(gather_if_gpu(layer_struct.Bias));
    else
        b = zeros(expected_size, 1);
    end
end
