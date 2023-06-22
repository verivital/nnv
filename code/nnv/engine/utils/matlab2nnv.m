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
    fprintf('\nParsing Layer %d... \n', i);
    customLayer_no_NLP = 0;
    try
       pat = digitsPattern(4);
       if exist(str2num(extract(string(L.Name),pat)),integer) && isempty(struct2array(L.ONNXParams.Nonlearnables))
        customLayer_no_NLP = 1;
       end
    catch
        try
            if contains(class(L), "PadLayer") && all(extractdata(struct2array(L.ONNXParams.Nonlearnables))==0)
                customLayer_no_NLP = 1;
            end
        catch
        end
    end
    

    % Layers with no effect on the reachability analysis
    if isa(L, 'nnet.cnn.layer.DropoutLayer') || isa(L, 'nnet.cnn.layer.SoftmaxLayer') || isa(L, 'nnet.cnn.layer.ClassificationOutputLayer') ...
            || isa(L,"nnet.onnx.layer.VerifyBatchSizeLayer") || isa(L, "nnet.cnn.layer.RegressionOutputLayer") || customLayer_no_NLP == 1
        Li = PlaceholderLayer.parse(L);

    % Image Input Layer
    elseif isa(L, 'nnet.cnn.layer.ImageInputLayer')
        Li = ImageInputLayer.parse(L);

    % Sequence Input Layer
    elseif isa(L, 'nnet.cnn.layer.SequenceInputLayer')
            Li = SequenceInputLayer.parse(L);

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

    % Max Pooling 2D Layer
    elseif isa(L, 'nnet.cnn.layer.MaxPooling2DLayer')
        Li = MaxPooling2DLayer.parse(L);

    % Average Pooling 2D Layer
    elseif isa(L, 'nnet.cnn.layer.AveragePooling2DLayer')
        Li = AveragePooling2DLayer.parse(L);
   
    % Global Average Pooling 2D Layer
    elseif isa(L, 'nnet.cnn.layer.GlobalAveragePooling2DLayer')
        Li = GlobalAveragePooling2DLayer.parse(L);

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
            || isa(L, 'nnet.onnx.layer.FlattenInto2dLayer')
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
        pairedMaxPoolingName = NN.getPairedMaxPoolingName(connects, Li.Name);
        Li.setPairedMaxPoolingName(pairedMaxPoolingName);

    % Depth Concatenation Layer (common in uNets)
    elseif isa(L, 'nnet.cnn.layer.DepthConcatenationLayer')
        Li = DepthConcatenationLayer.parse(L);
    
    % Reshape Layer (custom created after parsing ONNX layers)
    elseif contains(class(L), "ReshapeLayer")
        Li = ReshapeLayer.parse(L);

    % Reshape Layer (custom created after parsing ONNX layers)
    elseif contains(class(L), "UpsampleLayer")
        Li = UpsampleLayer.parse(L);
    
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

% Assigning layer names to correspnding index
name2number = containers.Map(names,indxs);

% ConnectionsTable = table(new_sources, new_dests, 'VariableNames', {'Source', 'Destination'});

% Create neural network
net = NN(nnvLayers, conns);
net.name2indx = name2number;

end