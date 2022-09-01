function controller = load_docking_controller(onnxfile)
    model = importONNXLayers(onnxfile);
    % 1) pre-processing (matmul 1)
    % 2) tanh layer (256x256?)
    % 3) post-processing (matmul too?)
    
    % Only "layer" with parameters is #2
    mlp = model.Layers(2);
    % How is the normalization playing out? Just multiply the inputs by those
    % values?
    pre = extractdata(mlp.model_preprocess_nor);
    post = extractdata(mlp.model_postprocess_fi); % This matrix won't work, need to add 2 rows of zeros to make it a NN function
    post(3:4,:) = zeros(2,4);
    
    % Weights and bias (looks like activation functions are linear, tanh,
    % linear)
    w1 = extractdata(mlp.model_model_fc_1_Mat);
    w2 = extractdata(mlp.model_model_fc_2_Mat);
    w3 = extractdata(mlp.model_model_fc_out_M);
    b1 = extractdata(mlp.ONNXParams.Nonlearnables.model_model_fc_1_Bia);
    b2 = extractdata(mlp.ONNXParams.Nonlearnables.model_model_fc_2_Bia);
    b3 = extractdata(mlp.ONNXParams.Nonlearnables.model_model_fc_out_B);
    
    %% Create Neural Network for NNV
    if exist('LayerS','class')
        % n = 5; % pre + 3 layers + post
        Lpre = LayerS(pre,zeros(4,1),'purelin');
        L1 = LayerS(w1,b1,'tansig');
        L2 = LayerS(w2,b2,'tansig');
        L3 = LayerS(w3,b3,'purelin');
        Lpost = LayerS(post,zeros(4,1),'tansig');
        Layers = [Lpre,L1,L2,L3,Lpost];
        controller = FFNNS(Layers);
    else
        controller = [];
        warning('NNV is not installed or added to the MATLAB path. Controller cannot be loaded, it will only be converted to NNV and ONNX is these files do not exist.')
    end

    %% Create file to match other NNV-saved NNs
    if ~isfile('model.mat')
        W = {}; b = {};
        % weights
        W{1} = pre; W{2} = w1; W{3} = w2; W{4} = w3; W{5} = post(1:2,:);
        % bias
        b{1} = zeros(4,1); b{2} = b1; b{3} = b2; b{4} = b3; b{5} = zeros(2,1);
        % activation functions
        act_fcns = ['linear';'tanh  ';'tanh  '; 'linear'; 'linear'];
        % Save model
        save('model.mat','W','b','act_fcns');
    end
    if ~isfile('bias_model.onnx')
        ToONNX('model.mat','bias_model.onnx');
    end
end