function [result, vTime] = run_vnncomp_instance(category, onnx, vnnlib, outputfile)
% single script to run all instances (supported) from the vnncomp2023

t = tic; % start timer
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
inputSize = net.Layers(1, 1).InputSize;
% disp(net);
% disp(nnvnet);

% Load property to verify
property = load_vnnlib(vnnlib);
lb = property.lb; % input lower bounds
ub = property.ub; % input upper bounds
prop = property.prop; % output spec to verify
% disp(property);
% disp(property.prop{1});

% Steps 2 and 3 may be executed differently depends on how the vnnlib
% properties are defined

%% 2) SAT?

nRand = 1000; % number of random inputs

% Choose how to falsify based on vnnlib file
if ~isa(lb, "cell") && length(prop) == 1 % one input, one output 
    counterEx = falsify_single(lb, ub, nRand, prop.Hg);
% elseif isa(lb, "cell") && length(lb) == length(prop) % multiple inputs, multiple outputs
%     
% elseif isa(lb, "cell") && length(prop) == 1 % can violate the output property from any of the input sets

else
    error("Working on adding support to other vnnlib properties")
end


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


%% 4) Process results

vTime = toc(t); % save total computation time

% Write results to output file
if result == 0
    fid = fopen(outputfile, 'w');
    fprintf(fid, 'sat \n');
    fclose(fid);
    write_counterexample(outputfile, counterEx)
elseif result == 1
    fid = fopen(outputfile, 'w');
    fprintf(fid, 'unsat \n');
    fclose(fid);
elseif result == 2
    fid = fopen(outputfile, 'w');
    fprintf(fid, 'unknown \n');
    fclose(fid);
end

quit; % does this work when running matlab.engine from python in background?


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

    elseif contains(category, "yolo")
        % yolo: onnx to matlab
        net = importONNXNetwork(onnx); % padlayer
        nnvnet = matlab2nnv(net);

    elseif contains(category, "acasxu")
        % acasxu: onnx to nnv
        net = importONNXNetwork(onnx, "InputDataFormats","BCSS");
        nnvnet = matlab2nnv(net);
        
    else % all other benchmarks
        % traffic: onnx to matlab: opset15 issues
        error("ONNX model not supported")
    end

end

function xRand = create_random_examples(lb, ub, nR, inputSize)
    xB = Box(lb, ub); % lb, ub must be vectors
    xRand = xB.sample(nR);
    xRand = reshape(xRand,[inputSize nR]); % reshape vectors into net input size
    xRand(:,:,:,nR+1) = x; % add lower bound 
    xRand(:,:,:,nR+2) = x; % add upper bound
end

function write_counterexample(outputfile, counterEx)
    % First line - > sat
    % after that, write the variables for each input dimension  of the counterexample
    %
    % Example:
    %  ( (X_0 0.12132)
    %    (X_1 3.45454)
    %    ( .... )
    %    (Y_0 2.32342)
    %    (Y_1 3.24355)
    %    ( ... )
    %    (Y_N 0.02456))
    %

end

function counterEx = falsify_single(lb, ub, nRand, Hs)
    counterEx = nan;
    xRand = create_random_examples(lb, ub, nRand, inputSize);
    s = size(xRand);
    n = length(z);
    %  look for counterexamples
    for i=1:s(n)
        x = get_example(xRand, i);
        yPred = predict(net, x);
        % check if property violated
        yPred = reshape(yPred, [], 1); % convert to column vector (if needed)
        for h=1:length(Hs)
            if Hs.contains(yPred) % property violated
                counterEx = {x; yPred}; % save input/output of countex-example
                break;
            end
        end
    end
    
end

function x = get_example(xRand,i)
    s = size(xRand);
    n = length(z);
    if n == 4
        x = xRand(:,:,:,i);
    elseif n == 3
        x = xRand(:,:,i);
    elseif n == 2
        x = xRand(:,i);
    else
        error("InputSize = "+string(s));
    end
end

