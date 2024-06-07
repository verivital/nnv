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
        elseif contains(class(L), "Reshape_To_ReshapeLayer") && length(fields(L.ONNXParams.Nonlearnables))==2
            % future: need to check if the previous layers output size == last reshape layers dimension
            customLayer_no_NLP = 1;
        end
    catch
        
    end
    

    % Layers with no effect on the reachability analysis
    if isa(L, 'nnet.cnn.layer.DropoutLayer') || isa(L, 'nnet.cnn.layer.SoftmaxLayer') || isa(L, 'nnet.cnn.layer.ClassificationOutputLayer') ...
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

    % Batch Normalization Layer
    elseif isa(L, 'nnet.cnn.layer.BatchNormalizationLayer')
        Li = BatchNormalizationLayer.parse(L);

    % Layer Normalization Layer
    elseif isa(L, 'nnet.cnn.layer.LayerNormalizationLayer')
        Li = LayerNormalizationLayer.parse(L);

    % Max Pooling 2D Layer
    elseif isa(L, 'nnet.cnn.layer.MaxPooling2DLayer')
        Li = MaxPooling2DLayer.parse(L);

    % Average Pooling 2D Layer
    elseif isa(L, 'nnet.cnn.layer.AveragePooling2DLayer')
        Li = AveragePooling2DLayer.parse(L);

    % Average Pooling 3D Layer
    elseif isa(L, 'nnet.cnn.layer.AveragePooling3DLayer')
        Li = AveragePooling3DLayer.parse(L);
   
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
    elseif isa(L, 'nnet.onnx.layer.ElementwiseAffineLayer')
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
    
    % Custom flatten layers (avoid if possible)
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

    % All other layers are currently not supported in NNV
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
