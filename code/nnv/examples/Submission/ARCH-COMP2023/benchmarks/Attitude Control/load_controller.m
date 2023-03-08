function controller = load_controller(onnxfile)
    model = importONNXLayers(onnxfile);
    
    % Only "layer" with parameters is #2
    mlp = model.Layers(2);
    
    % Weights and bias (looks like activation functions are linear, tanh,
    % linear)
    w1 = extractdata(mlp.layers_lin1_weight)';
    w2 = extractdata(mlp.layers_lin2_weight)';
    w3 = extractdata(mlp.layers_lin3_weight)';
    w4 = extractdata(mlp.layers_lin4_weight)';
    b1 = extractdata(mlp.layers_lin1_bias);
    b2 = extractdata(mlp.layers_lin2_bias);
    b3 = extractdata(mlp.layers_lin3_bias);
    b4 = extractdata(mlp.layers_lin4_bias);
    
    %% Create Neural Network for NNV
    if exist('LayerS','class')
        % n = 4; % 3 hidden sigmoid layers + linear output
        L1 = LayerS(w1,b1,'logisg');
        L2 = LayerS(w2,b2,'logsig');
        L3 = LayerS(w3,b3,'logsig');
        L4 = LayerS(w4,b4,'purelin');
        Layers = [L1,L2,L3,L4];
        controller = FFNNS(Layers);
    else
        controller = [];
        warning('NNV is not installed or added to the MATLAB path. Controller cannot be loaded, it will only be converted to NNV and ONNX is these files do not exist.')
    end

    %% Create file to match other NNV-saved NNs
    if ~isfile('model.mat')
        W = {}; b = {};
        % weights
        W{1} = w1; W{2} = w2; W{3} = w3; W{4} = w4;
        % bias
        b{1} = b1; b{2} = b2; b{3} = b3; b{4} = b4;
        % activation functions
        act_fcns = ['logsig';'logsig';'logsig'; 'linear'];
        % Save model
        save('model.mat','W','b','act_fcns');
    end
    if ~isfile('model.onnx')
        ToONNX('model.mat','model.onnx');
    end
end