function nn = load_NN_from_mat(matfile)
% nn = loadNN_from_mat(matfile) 
% 
% - Output:
%     nn = NN object (common class for all deep learning models in NNV)
% - Input
%     matfile = format used in old NNV, which has 3 variables, W -> weights, 
%              b -> bias, act_fcns -> activation function for each layer
% This format was only used for fullyconnected networks (FFNNS), now deprecated

    % load data
    net_info = load(matfile);
    % Check all data is there
    if ~isfield(net_info,'W') || ~isfield(net_info, 'b') || ~isfield(net_info, 'act_fcns')
        error("Wrong matfile format, W, b or act_fcns is missing")
    end

    % Extract neural network data
    W = net_info.W;          % weights
    b = net_info.b;          % bias
    n = length(b);           % number of (weighted) layers
    acf = net_info.act_fcns; % activation functions

    % Check for syntax correctness 
    if ~isa(W, "cell") || ~isa(b, "cell") || ~isa(acf, "char")
        error("Wrong data type on one or more of the network variables")
    end
    
    % Create array of NN layers
    Layers = {}; % cell array to store layers
    layer_count = 1;
    for i=1:n
        Layers{layer_count} = FullyConnectedLayer(W{i},b{i});
        layer_count = layer_count + 1;
        layer  = ActFunction(acf(i,:));
        if ~isstring(layer) % if not linear or identity, then add layer
            Layers{layer_count} = layer;
            layer_count = layer_count + 1;
        end
    end
    
    % Get output size of network from last fullyconnected layer
    for l = length(Layers):-1:1
        if isa(Layers{l}, "FullyConnectedLayer")
            outputSize = Layers{l}.OutputSize;
        end
    end

    % Create the NN object
    nn = NN(Layers, Layers{1}.InputSize, outputSize, 'nn'); % cell array of layer, connection table, input size, output size, name
end

%% Helper function
function layer = ActFunction(act)
% Generate NN layers based of activation functions

    act = strtrim(act); % Remove empty white spaces
    if contains(act,'relu1')
        layer = SaturatingLinearLayer();
    elseif contains(act,'linear')
        layer = "identity"; % no need to add any layer
    elseif contains(act,'relu2')
        layer = SaturatingLinearSymmLayer();
    elseif contains(act,'relu')
        layer = ReluLayer();
    elseif contains(act,'tanh') || contains(act, 'tansig')
        layer = TanhLayer();
    elseif contains(act,'igmoid') || contains(act, 'logsig')
        layer = SigmoidLayer;
    else
        error("Unkown activation function ("+act+")");
    end
    
end
