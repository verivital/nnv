function net = matlab2nnv(Mnetwork)

%% Check for correct inputs and process network
ntype = class(Mnetwork); % input type
if ~contains(ntype, ["SeriesNetwork", "LayerGraph", "DAGNetwork", "dlnetwork"])
    error('Wrong input type. Input must be a SeriesNetwork, LayerGraph, DAGNetwork, or dlnetwork');
end

%% Process input types
% Input is a MATLAB type neural network (layergraph, seriesNetwork, dlnetwork or dagnetwork)
if ntype== "SeriesNetwork"
    conns = layerGraph(Mnetwork).Connections; % get the table of connections
else
    conns = Mnetwork.Connections; % get the table of connections
end

Layers = Mnetwork.Layers; % get the list of layers

%% Transform to NNV 

n = length(Layers);
nnvLayers = cell(n,1);
names = strings(n,1);

% Parse network layer-by-layer
for i=1:n
    L = Layers(i);
%     fprintf('\nParsing Layer %d... \n', i);
    customLayer_no_NLP = 0;
    try
        if contains(class(L), "FlattenLayer") && isempty(L.ONNXParams.Nonlearnables)
            customLayer_no_NLP = 1;
        elseif contains(class(L), "PadLayer") %&& all(extractdata(struct2array(L.ONNXParams.Nonlearnables))==0)
            customLayer_no_NLP = check_layer_parameters(L);
        elseif contains(class(L), "Reshape_To_ReshapeLayer")
            % A fused pure-Reshape chain (the class-name prefix encodes the fused
            % ops: Reshape_To_ReshapeLayer = Reshape->...->Reshape ONLY, no weights
            % -- unlike e.g. Gemm_To_ReshapeLayer). A double reshape whose FINAL
            % target is flat is a data no-op for the flat star pipeline.
            % future: check previous layer's output size == last reshape target dim
            if isprop(L, 'ONNXParams') && length(fields(L.ONNXParams.Nonlearnables))==2
                % pre-R2026a custom-layer form
                customLayer_no_NLP = 1;
            elseif isprop(L, 'Vars')
                % R2026a importNetworkFromONNX custom layers carry the ONNX
                % initializers in .Vars (no ONNXParams property). Accept the fused
                % double-reshape as a no-op ONLY when the final reshape target is
                % flat ([-1 N]): reshape-to-image then flatten-back preserves the
                % element order, so the layer is the identity on flat data. Any
                % other shape falls through to ReshapeLayer.parse (which fails
                % loud -> the runner reports unknown; sound).
                vf = fieldnames(L.Vars);
                if numel(vf) == 2
                    lastShape = L.Vars.(vf{end});
                    try, lastShape = extractdata(lastShape); catch, end
                    lastShape = double(lastShape(:)');
                    if numel(lastShape) == 2 && lastShape(1) == -1
                        customLayer_no_NLP = 1;
                    end
                end
            end
        end
    catch
        
    end
    

    % SOUNDNESS: a MID-NETWORK softmax must NOT become an identity placeholder
    % (identity is only sound for a FINAL softmax under argmax-on-logits specs).
    % "Final" = every subsequent layer has no effect on reachability; erring
    % toward non-final is sound (the [0,1] bounds are merely looser).
    if isa(L, 'nnet.cnn.layer.SoftmaxLayer') && ~is_final_softmax(Layers, i, conns)
        Li = SoftmaxLayer(L.Name);
        Li.IsFinalLayer = false;   % sound interval bounds (see SoftmaxLayer.reach)

    elseif contains(class(L), "LogSoftmax") && ~is_final_softmax(Layers, i, conns)
        % Log-softmax outputs lie in (-inf, 0] -- neither identity nor the
        % softmax [0,1] bounds is sound for it mid-network. Refuse loudly.
        error('matlab2nnv:midNetworkLogSoftmax', ...
            ['Mid-network LogSoftmax layer ''%s'' is not supported soundly ' ...
             '(identity passthrough would be unsound).'], L.Name);

    % Layers with no effect on the reachability analysis
    elseif isa(L, 'nnet.cnn.layer.DropoutLayer') || isa(L, 'nnet.cnn.layer.SoftmaxLayer') || isa(L, 'nnet.cnn.layer.ClassificationOutputLayer') ...
            || isa(L,"nnet.onnx.layer.VerifyBatchSizeLayer") || isa(L, "nnet.cnn.layer.RegressionOutputLayer") ...
            || customLayer_no_NLP == 1 || isa(L, "nnet.onnx.layer.CustomOutputLayer") || contains(class(L), "LogSoftmax")
        Li = PlaceholderLayer.parse(L);

    elseif customLayer_no_NLP
        Li = PlaceholderLayer.parse(L);

    % Image Input Layer
    elseif isa(L, 'nnet.cnn.layer.ImageInputLayer')
        Li = ImageInputLayer.parse(L);
    
    % Image 3D Input Layer 
    elseif isa(L, 'nnet.cnn.layer.Image3DInputLayer')
        Li = Image3DInputLayer.parse(L);

    % Sequence Input Layer
    elseif isa(L, 'nnet.cnn.layer.SequenceInputLayer')
            Li = SequenceInputLayer.parse(L);

    % Convolutional 3D layer
    elseif isa(L, 'nnet.cnn.layer.Convolution3DLayer') 
        Li = Conv3DLayer.parse(L);

    % Convolutional 2D layer
    elseif isa(L, 'nnet.cnn.layer.Convolution2DLayer') 
        Li = Conv2DLayer.parse(L);
    
    % Convolutional 1D layer
    elseif isa(L, 'nnet.cnn.layer.Convolution1DLayer') 
            Li = Conv1DLayer.parse(L);

     % Transposed Convolution 1D layer
    elseif isa(L, 'nnet.cnn.layer.TransposedConvolution1DLayer') 
            Li = TransposedConv1DLayer.parse(L);

    % ReLU Layer (also referred to as poslin)
    elseif isa(L, 'nnet.cnn.layer.ReLULayer')
        Li = ReluLayer.parse(L);

    % LeakyReLU Layer
    elseif isa(L, 'nnet.cnn.layer.LeakyReLULayer')
            Li = LeakyReluLayer.parse(L);

    % Tanh Layer
    elseif isa(L, 'nnet.cnn.layer.TanhLayer')
            Li = TanhLayer.parse(L);

    % GELU Layer (Gaussian Error Linear Unit - used in Transformers)
    elseif isa(L, 'nnet.cnn.layer.GELULayer')
            Li = GeluLayer.parse(L);

    % Batch Normalization Layer
    elseif isa(L, 'nnet.cnn.layer.BatchNormalizationLayer')
        Li = BatchNormalizationLayer.parse(L);

    % Layer Normalization Layer
    elseif isa(L, 'nnet.cnn.layer.LayerNormalizationLayer')
        Li = LayerNormalizationLayer.parse(L);

    % Self-Attention Layer (Transformer attention)
    elseif isa(L, 'nnet.cnn.layer.SelfAttentionLayer')
        Li = MultiHeadAttentionLayer.parse(L);

    % Max Pooling 2D Layer
    elseif isa(L, 'nnet.cnn.layer.MaxPooling2DLayer')
        Li = MaxPooling2DLayer.parse(L);

    % Average Pooling 2D Layer
    elseif isa(L, 'nnet.cnn.layer.AveragePooling2DLayer')
        Li = AveragePooling2DLayer.parse(L);

    % Average Pooling 3D Layer
    elseif isa(L, 'nnet.cnn.layer.AveragePooling3DLayer')
        Li = AveragePooling3DLayer.parse(L);
   
    % Global Average Pooling 3D Layer
    elseif isa(L, 'nnet.cnn.layer.GlobalAveragePooling3DLayer')
        Li = GlobalAveragePooling3DLayer.parse(L);

    % Global Average Pooling 2D Layer
    elseif isa(L, 'nnet.cnn.layer.GlobalAveragePooling2DLayer')
        Li = GlobalAveragePooling2DLayer.parse(L);

    % Global Average Pooling 1D Layer
    elseif isa(L, 'nnet.cnn.layer.GlobalAveragePooling1DLayer')
        Li = GlobalAveragePooling1DLayer.parse(L);

    % Fully Connected Layer
    elseif isa(L, 'nnet.cnn.layer.FullyConnectedLayer')
        Li = FullyConnectedLayer.parse(L);

    % Addition Layer
    elseif isa(L, 'nnet.cnn.layer.AdditionLayer')
        Li = AdditionLayer.parse(L);
    
    % Concatenation Layer
    elseif isa(L, 'nnet.cnn.layer.ConcatenationLayer')
        Li = ConcatenationLayer.parse(L);

    % Pixel Classification Layer (used for Semantic Segmentation output)
    elseif isa(L, 'nnet.cnn.layer.PixelClassificationLayer')
        Li = PixelClassificationLayer.parse(L);

    % Flatten Layer
    elseif isa(L, 'nnet.keras.layer.FlattenCStyleLayer') || isa(L, 'nnet.cnn.layer.FlattenLayer') || isa(L, 'nnet.onnx.layer.FlattenLayer') ...
            || isa(L, 'nnet.onnx.layer.FlattenInto2dLayer') || isa(L, 'nnet.onnx.layer.Flatten3dInto2dLayer')
        Li = FlattenLayer.parse(L);

    % Sigmoid Layer (also referred to as logsig)
    elseif isa(L, 'nnet.keras.layer.SigmoidLayer') || isa(L, 'nnet.onnx.layer.SigmoidLayer') || isa(L, 'nnet.cnn.layer.SigmoidLayer')
        Li = SigmoidLayer.parse(L);

    % ElementWise Affine Layer (often used as a bias layer after FC layers)
    % R2024b and earlier ONNX imports emit nnet.onnx.layer.ElementwiseAffineLayer.
    % R2025a+ ONNX imports emit nnet.cnn.layer.ScalingLayer with identical
    % Y = Scale.*X + Offset semantics (MathWorks release note:
    % "ONNX Add, Mul, Sub, Neg, Div imported as scalingLayer"). Both route to
    % the same parse() -- ElementwiseAffineLayer.parse now accepts either class.
    elseif isa(L, 'nnet.onnx.layer.ElementwiseAffineLayer') || isa(L, 'nnet.cnn.layer.ScalingLayer')
        % Both R2024b- ElementwiseAffineLayer and R2025a+ ScalingLayer route to
        % the single parse() (Y = Scale.*X + Offset). [27]: the previous code
        % had a SECOND, unreachable ScalingLayer branch below this one (this
        % branch already matches nnet.cnn.layer.ScalingLayer) whose hand-rolled
        % logic diverged from parse() -- forcing DoScale/DoOffset both true and
        % substituting a zero offset. Removed to keep one code path.
        Li = ElementwiseAffineLayer.parse(L);

    % Feature input layer
    elseif isa(L, 'nnet.cnn.layer.FeatureInputLayer') || isa(L, 'nnet.onnx.layer.FeatureInputLayer')
        Li = FeatureInputLayer.parse(L); 

    % Transposed Convolution 2D Layer
    elseif isa(L, 'nnet.cnn.layer.TransposedConvolution2DLayer')
        Li = TransposedConv2DLayer.parse(L);

    % Max Unpooling 2D Layer
    elseif isa(L, 'nnet.cnn.layer.MaxUnpooling2DLayer')
        Li = MaxUnpooling2DLayer.parse(L, conns);

    % Depth Concatenation Layer (common in uNets)
    elseif isa(L, 'nnet.cnn.layer.DepthConcatenationLayer')
        Li = DepthConcatenationLayer.parse(L);
    
    % Reshape Layer (custom created after parsing ONNX layers)
    elseif contains(class(L), "ReshapeLayer")
        Li = ReshapeLayer.parse(L);

    % Reshape Layer (custom created after parsing ONNX layers)
    elseif contains(class(L), "UpsampleLayer")
        Li = UpsampleLayer.parse(L);

    % Reshape Layer (custom created after parsing ONNX layers)
    elseif isa(L, 'nnet.cnn.layer.LSTMLayer')
        Li = LstmLayer.parse(L);

    elseif isa(L, 'nnet.cnn.layer.Resize2DLayer')
        Li = Resize2DLayer.parse(L);

    elseif contains(class(L), 'Reshape_To_ConcatLayer')
        Li = ReshapeToConcatenationLayer.parse(L);
    
    elseif contains(class(L), 'MatMul_To_AddLayer')
        Li = MatMulToAddLayer.parse(L);

    elseif contains(class(L), 'MatMul_To_ReduceSumLayer')
        Li = MatMulToReduceSumLayer.parse(L);

    elseif contains(class(L), 'MatMul_To_SubLayer')
        Li = MatMulToSubLayer.parse(L);

    elseif contains(class(L), ["flatten"; "Flatten"])
        % Check previous layer to see if we can neglect this one in the analysis
        for k=i-1:-1:1
            if contains(class(nnvLayers{k}), 'Input')
                if ~strcmp(nnvLayers{k}.Normalization, 'none')
                    fprintf('Layer %d is a %s which have not supported yet in nnv, please consider removing this layer for the analysis \n', i, class(L));
                    error('Unsupported Class of Layer');
                end
            elseif ~isa(nnvLayers{k}, 'PlaceholderLayer')
                fprintf('Layer %d is a %s which have not supported yet in nnv, please consider removing this layer for the analysis \n', i, class(L));
                error('Unsupported Class of Layer');
            end
        end
        % If we can neglect all previous layers, reinitialize layers and parse them again as placeholder layers
        nnvLayers = cell(n,1);
        % Parse all previous layers again
        for li = 1:i-1
            L = Layers(li);
            Li = PlaceholderLayer.parse(L);
            nnvLayers{li} = Li;
        end
        % Parse current flatten layer
        L = Layers(i);
        Li = PlaceholderLayer.parse(L);

        % All other layers are currently not supported in NN
    else
        fprintf('Layer %d is a %s which have not supported yet in nnv, please consider removing this layer for the analysis \n', i, class(L));
        error('Unsupported Class of Layer');
    end

    % Add layer name to list
    names(i) = string(L.Name);
    nnvLayers{i} = Li;
end
indxs = 1:n;


% We compute the reachability and evaluation layer by layer, executing the
% connections in the order they appear, so need to ensure the order makes sense:
%  - sometimes some of the skipped connections (like in resnet or unets), they appear at the end so NNV returns the wrong output

% Revert back to previous working code. Navid's changes were getting stuck
% on the while loop...

if height(conns) == length(nnvLayers) - 1 % fullyconnected layers, no skips
    % Assigning layer names to correspnding index
    name2idx = containers.Map(names,indxs);
    nnvConns = conns;
else
    [nnvLayers, nnvConns, name2idx] = process_connections(nnvLayers, conns, names, indxs);
end

% ConnectionsTable = table(new_sources, new_dests, 'VariableNames', {'Source', 'Destination'});

% Create neural network
net = NN(nnvLayers, nnvConns);
net.name2indx = name2idx;

% Get input and output sizes
outSize = getOutputSize(Mnetwork);
net.OutputSize = outSize; 

net.matlabnet = Mnetwork;
if ntype== "SeriesNetwork"
    net.matlabnet = dag2dlnetwork(net.matlabnet);
end

end


% Helper function for connections
function [nnvLayers, nnvConns, name2number] = process_connections(nnvLayers, conns, names, idxs)
    % Assigning layer names to correspnding index
    name2number = containers.Map(names,idxs);

    % Step 1 - initialize a variable to keep track of the number of inputs in NNV
    count_inputs = ones(length(nnvLayers),1); % order is same as nnvLayers (default = 1)
    for k=1:length(nnvLayers)
        if contains(class(nnvLayers{k}), ["Concatenation", "Addition"])
            count_inputs(k) = nnvLayers{k}.NumInputs;
        elseif contains(class(nnvLayers{k}), "Input")
            count_inputs(k) = 0; % this is the initial layer, no prior connections
        end
    end

    % Step 2 - ensure these number of inputs are seen before executing
    %       that layer, otherwise, move to the next connection
    final_sources = [];
    final_dests = [];
    orig_sources = conns.Source;
    orig_dests = conns.Destination;
    skipped = 1: height(conns);%idxs; %start with assumption that all connections are skipped

    % start adding cnnections to the final list
    while ~isempty(skipped) 
        % Get connection idx
        c = skipped(1);
        skipped(1) = []; % like popping the first connection from the list
        % get connection names
        source_name = orig_sources{c};
        source_layer_name = split(source_name, '/');
        if length(source_layer_name) > 1 % either a layer with multiple dets or some "properties" like size are getting sent to other layers
            sec_name = source_layer_name{2};
            if ~contains(sec_name, 'out')
                continue;
            end
        end
        source_layer_name = source_layer_name{1};
        dest_name = orig_dests{c};
        % parse dest_name as it could be defined as name/in1 and so on, then add it back at the end
        dest_layer_name = split(dest_name, '/');
        dest_layer_name = dest_layer_name{1};

        % check if source has an input counter of 0
        source_idx = name2number(source_layer_name);
        if count_inputs(source_idx) == 0
            dest_idx = name2number(dest_layer_name);
            count_inputs(dest_idx) = count_inputs(dest_idx)-1;
            final_sources = [final_sources; string(source_layer_name)];
            final_dests = [final_dests; string(dest_name)];
        else
            if isempty(skipped)
                error("There is an unreachable connection")
            else
                skipped = [skipped, c];
            end
        end        
        
    end

    nnvConns = table(final_sources, final_dests, 'VariableNames', {'Source', 'Destination'});    

end


% Check if custom layers can be skipped (e.g., pad layer with 0 padding does not change anything in the network)
function status = check_layer_parameters(L)
    status = 0;
    if isempty(L.ONNXParams.Nonlearnables)
        status = 1;
    elseif contains(class(L), 'PadLayer')
        params = struct2cell(L.ONNXParams.Nonlearnables);
        for i=1:length(params)
            tmp = extractdata(params{i});
            if ~all(tmp==0)
                error('Layer not supported');
            end
        end
        status = 1;
    end
end

% A softmax is "final" (identity placeholder is sound for argmax-on-logits
% specs) only if every GRAPH SUCCESSOR has no effect on reachability. Erring
% toward non-final is SOUND (bounds get looser, never wrong); only treating a
% genuinely mid-network softmax as final would be unsound.
%
% [22]: decide by Connections topology, not array order. The Layers array is not
% guaranteed to be topologically sorted, so walking Layers(idx+1:end) by array
% position could call a softmax 'final' when a real layer follows it in the
% GRAPH but precedes it in the array. When a usable Connections table is
% available we BFS the actual successors; otherwise we fall back to the array
% walk (sound for the usual topologically-ordered arrays).
function tf = is_final_softmax(Layers, idx, conns)
    if nargin >= 3 && istable(conns) && ~isempty(conns) && ...
            all(ismember({'Source','Destination'}, conns.Properties.VariableNames))
        tf = is_final_softmax_topo(Layers, idx, conns);
        return;
    end
    tf = true;
    for k = idx+1:numel(Layers)
        if ~softmax_successor_inert(Layers(k))
            tf = false; return;
        end
    end
end

% True iff a successor layer has no effect on reachability (so a softmax feeding
% only such layers is "final").
function tf = softmax_successor_inert(Lk)
    tf = isa(Lk, 'nnet.cnn.layer.ClassificationOutputLayer') || ...
         isa(Lk, 'nnet.cnn.layer.RegressionOutputLayer') || ...
         isa(Lk, 'nnet.cnn.layer.DropoutLayer') || ...
         isa(Lk, "nnet.onnx.layer.CustomOutputLayer") || ...
         isa(Lk, "nnet.onnx.layer.VerifyBatchSizeLayer");
end

% Topology-aware finality: BFS the softmax's graph successors via Connections.
% Final iff EVERY reachable successor layer is inert. Layer names in the conns
% table may carry '/port' suffixes -- strip them to match Layers(k).Name.
function tf = is_final_softmax_topo(Layers, idx, conns)
    names = strings(numel(Layers), 1);
    for k = 1:numel(Layers)
        names(k) = string(Layers(k).Name);
    end
    src = strip_port(string(conns.Source));
    dst = strip_port(string(conns.Destination));
    start = string(Layers(idx).Name);

    tf = true;
    visited = false(numel(Layers), 1);
    frontier = start;
    while ~isempty(frontier)
        nextNodes = strings(0,1);
        for f = 1:numel(frontier)
            succ = dst(src == frontier(f));   % direct graph successors
            for s = 1:numel(succ)
                li = find(names == succ(s), 1);
                if isempty(li) || visited(li), continue; end
                visited(li) = true;
                if ~softmax_successor_inert(Layers(li))
                    tf = false; return;       % a real layer follows -> not final
                end
                nextNodes(end+1,1) = succ(s); %#ok<AGROW>
            end
        end
        frontier = nextNodes;
    end
end

function s = strip_port(s)
    % 'layerName/out1' -> 'layerName'
    s = extractBefore(s + "/", "/");
end

% Get output size of the network
function outSize = getOutputSize(network)
    outSize = 0; % default
    if isa(network, "SeriesNetwork")
        outputLayer = network.Layers(end);
        if isa(outputLayer, "nnet.cnn.layer.ClassificationOutputLayer")
            outSize = outputLayer.OutputSize;
        end
    end
end
