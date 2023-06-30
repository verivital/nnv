function result = run_vnncomp_instance(category, onnx, vnnlib)
% single script to run all instances (supported) from the vnncomp2023

result = 2; % unknown (to start with)

% Process:
%  1) Load components
%  2) SAT? - Randomly evaluate network to search for counterexample 
%  3) UNSAT? 
%     a) Compute reachability 
%     b) Verify (compute intersection between reach set and output property halfspace)

%% 1) Load components

% Load networks
[net, nnvnet] = load_vnncomp_network(category, onnx);
disp(net);
disp(nnvnet);

% Load property to verify
[lb,ub,prop] = load_vnnlib(vnnlib);
disp(prop);

% Steps 2 and 3 may be executed differently depends on how the vnnlib
% properties are defined

%% 2) SAT?

% Choose how to falsify based on vnnlib file
% if ~isa(lb, "cell") && length(prop) == 1 % one input, one output 
% 
% elseif isa(lb, "cell") && length(lb) == length(prop) % multiple inputs, multiple outputs
% 
% elseif isa(lb, "cell") && length(prop) == 1
% 
% else
%     error("Working on adding support to other vnnlib properties")
% end


%% 3) UNSAT?

% Choose how to verify based on vnnlib file
% if ~isa(lb, "cell") && length(prop) == 1 % one input, one output 
% 
% elseif isa(lb, "cell") && length(lb) == length(prop) % multiple inputs, multiple outputs
% 
% elseif isa(lb, "cell") && length(prop) == 1
% 
% else
%     error("Working on adding support to other vnnlib properties")
% end

end


function [net,nnvnet] = load_vnncomp_network(category, onnx)
% load vnncomp 2023 benchmark NNs (subset support)

    % collins_rul: onnx to nnvnet
    % collins_nets = load_collins_NNs;
    if contains(category, 'collins_rul')
        net = importONNXNetwork(onnx);
        nnvnet = matlab2nnv(onnx);

    elseif contains(category, "nn4sys")
        % nn4sys: onnx to matlab:
        net = importONNXLayers(onnx, "OutputDataFormats", "BC"); % lindex
        nnvnet = matlab2nnv(net);
        
    elseif contains(category, "dist_shift")
        % dist_shift: onnx to matlab:
        net = importONNXNetwork(onnx, "InputDataFormats", "BC", 'OutputDataFormats',"BC"); %reshape
        nnvnet = "";
        
    elseif contains(category, "cgan")
        % cgan
        net = importONNXNetwork(onnx,"InputDataFormats", "BC", 'OutputDataFormats',"BC"); %reshape
        nnvnet = matlab2nnv(net);
        
    elseif contains(category, "vgg16")
        % vgg16: onnx to matlab
        net = importONNXNetwork('vgg16-7.onnx'); % flattenlayer
        %reshapedInput = python_reshape(input,net_vgg.Layers(1,1).InputSize); % what is the input? assume it's all the same?
        %nnvnet = matlab2nnv(net);
        nnvnet = "";
        
    elseif contains(category, "tllverify")
        % tllverify: onnx to matlab
        net = importONNXNetwork(onnx,"InputDataFormats", "BC", 'OutputDataFormats',"BC");
        nnvnet = matlab2nnv(net);
        
    elseif contains(category, "vit")
        % vit: onnx to matlab
        net = importONNXNetwork(onnx, "TargetNetwork","dlnetwork" );
        nnvnet = "";
        
    elseif contains(category, "cctsdb_yolo")
        % cctsdb_yolo: onnx to matnet
        net = importONNXNetwork(onnx, "TargetNetwork","dlnetwork" );
        nnvnet = "";
        
    elseif contains(category, "collins_yolo")
        % collins_yolo: onnx to matlab:
        net = importONNXNetwork(onnx, "TargetNetwork","dlnetwork" );
        nnvnet = "";

    elseif contains(category, yolo)
        % yolo: onnx to matlab
        net = importONNXNetwork(onnx); % padlayer
        nnvnet = matlab2nnv(net);
        
    else % all other benchmarks
        % traffic: onnx to matlab: opset15 issues
        error("ONNX model not supported")
    end

end