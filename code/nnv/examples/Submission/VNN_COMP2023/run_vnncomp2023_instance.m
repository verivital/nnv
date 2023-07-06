function [status, vTime] = run_vnncomp2023_instance(category, onnx, vnnlib, outputfile)
% single script to run all instances (supported) from the vnncomp2023

t = tic; % start timer
status = 2; % unknown (to start with)

% Process:
%  1) Load components
%  2) SAT? - Randomly evaluate network to search for counterexample 
%  3) UNSAT? 
%     a) Compute reachability 
%     b) Verify (compute intersection between reach set and output property halfspace)

%% 1) Load components

% Load networks
% have to go to this path for the networks to load properly... lovely
old_path = pwd;
% cd /home/dieman95/Documents/MATLAB/nnv/code/nnv/examples/Submission/VNN_COMP2023/networks2023/;
cd /home/ubuntu/toolkit/code/nnv/examples/Submission/VNN_COMP2023/networks2023/;
[net, nnvnet] = load_vnncomp_network(category, onnx);
inputSize = net.Layers(1, 1).InputSize;
% disp(net);
% disp(nnvnet);

% Load property to verify
warning(vnnlib); % somehow it is failing to open the vnnlib files...
property = load_vnnlib(vnnlib);
lb = property.lb; % input lower bounds
ub = property.ub; % input upper bounds
prop = property.prop; % output spec to verify
% disp(property);
% disp(property.prop{1});

% Steps 2 and 3 may be executed differently depends on how the vnnlib
% properties are defined

%% 2) SAT?

nRand = 10; % number of random inputs

% Choose how to falsify based on vnnlib file
if ~isa(lb, "cell") && length(prop) == 1 % one input, one output 
    counterEx = falsify_single(net, lb, ub, inputSize, nRand, prop{1}.Hg);
% elseif isa(lb, "cell") && length(lb) == length(prop) % multiple inputs, multiple outputs
%     
% elseif isa(lb, "cell") && length(prop) == 1 % can violate the output property from any of the input sets

else
    warning("Working on adding support to other vnnlib properties")
end

% disp(toc(t));


%% 3) UNSAT?

% Define reachability options
reachOptions = struct;
% reachOptions.reachMethod = 'exact-star';
% reachOptions.reachOption = 'parallel';
% reachOptions.numCores = feature('numcores');
reachOptions.reachMethod = 'approx-star';

if status == 2 && isa(nnvnet, "NN")  && ~isa(reachOptions, 'struct') % no counterexample found and supported for reachability (otherwise, skip step 3 and write results)

% Choose how to verify based on vnnlib file
    if ~isa(lb, "cell") && length(prop) == 1 % one input, one output 
        % Get input set
        lb = reshape(lb, inputSize);
        ub = reshape(ub, inputSize);
        IS = ImageStar(lb, ub);
        % Compute reachability
        ySet = nnvnet.reach(IS, reachOptions);
        status = verify_specification(ySet, prop);
%     elseif isa(lb, "cell") && length(lb) == length(prop) % multiple inputs, multiple outputs
%     
%     elseif isa(lb, "cell") && length(prop) == 1
    
    else
        warning("Working on adding support to other vnnlib properties")
    end

end

%% 4) Process results

vTime = toc(t); % save total computation time

disp(status);
disp(vTime);
disp( " ");

% Write results to output file
if status == 0
    fid = fopen(outputfile, 'w');
    fprintf(fid, 'sat \n');
    fclose(fid);
    write_counterexample(outputfile, counterEx)
elseif status == 1
    fid = fopen(outputfile, 'w');
    fprintf(fid, 'unsat \n');
    fclose(fid);
elseif status == 2
    fid = fopen(outputfile, 'w');
    fprintf(fid, 'unknown \n');
    fclose(fid);
end

cd(old_path); % go back to where we ran the functions from
% quit; % does this work when running matlab.engine from python in background?


end

function [net,nnvnet] = load_vnncomp_network(category, onnxFile)
% load vnncomp 2023 benchmark NNs (subset support)

    onnx = split(onnxFile,'/');
    if iscell(onnx)
        onnx = onnx{end}; % cell
    else
        onnx = char(onnx(end)); % string
    end
    onnx = [onnx(1:end-5), '.mat'];

    % collins_rul: onnx to nnvnet
    % collins_nets = load_collins_NNs;
    if contains(category, 'collins_rul')
        net = load(onnx);
        net = net.net;
        nnvnet = matlab2nnv(net);

    elseif contains(category, "ml4acopf")
        net = load(onnx);
        net = net.net;
        nnvnet = "";

    elseif contains(category, "nn4sys")
        % nn4sys: onnx to matlab
        if ~contains(onnxFile, "2048")
            net = load(onnx);
            net = net.net;
            nnvnet = matlab2nnv(net);
        else
            error("We don't have those");
        end
        
    elseif contains(category, "dist_shift")
        % dist_shift: onnx to matlab:
        net = load(onnx);
        net = net.net;        
        nnvnet = "";
        
    elseif contains(category, "cgan")
        % cgan
        if ~contains(onnxFile, 'transformer')
            net = load(onnx);
            net = net.net;
            nnvnet = matlab2nnv(net);
        else
            error("Loading this one is not supported...")
        end
        
    elseif contains(category, "vgg")
        % vgg16: onnx to matlab
        net = load(onnx);
        net = net.net;
        %reshapedInput = python_reshape(input,net_vgg.Layers(1,1).InputSize); % what is the input? assume it's all the same?
        %nnvnet = matlab2nnv(net);
        nnvnet = "";
        
    elseif contains(category, "tllverify")
        % tllverify: onnx to matlab
        net = load(onnx);
        net = net.net;
        nnvnet = matlab2nnv(net);
        
    elseif contains(category, "vit")
        % vit: onnx to matlab
        net = load(onnx);
        net = net.net;
        nnvnet = "";
        
    elseif contains(category, "cctsdb_yolo")
        % cctsdb_yolo: onnx to matnet
        net = load(onnx);
        net = net.net;
        nnvnet = "";
        
    elseif contains(category, "collins_yolo")
        % collins_yolo: onnx to matlab:
        net = load(onnx);
        net = net.net;
        nnvnet = "";

    elseif contains(category, "yolo")
        % yolo: onnx to matlab
        net = load(onnx);
        net = net.net;
        nnvnet = matlab2nnv(net);

    elseif contains(category, "acasxu")
        % acasxu: onnx to nnv
        net = load(onnx);
        net = net.net;
        nnvnet = matlab2nnv(net);

    elseif contains(category, 'test')
        % test: onnx to nnv
        if contains(onnx, "sat")
            net = load(onnx);
            net = net.net;
            nnvnet = matlab2nnv(net);
        else
            error("Loading this one is not supported...")
        end

    else % all other benchmarks
        % traffic: onnx to matlab: opset15 issues
        error("ONNX model not supported")
    end

end

% Had to change this because importer fails using the command line...s
function [net,nnvnet] = load_vnncomp_network_local(category, onnx)
% load vnncomp 2023 benchmark NNs (subset support)

    % collins_rul: onnx to nnvnet
    % collins_nets = load_collins_NNs;
    if contains(category, 'collins_rul')
        net = importONNXNetwork(onnx);
        nnvnet = matlab2nnv(net);

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

function xRand = create_random_examples(net, lb, ub, nR, inputSize)
    xB = Box(lb, ub); % lb, ub must be vectors
    xRand = xB.sample(nR-2);
    xRand = [lb, ub, xRand];
    xRand = reshape(xRand,[inputSize nR]); % reshape vectors into net input size
    if isa(net, 'dlnetwork') % need to convert to dlarray
        if isa(net.Layers(1, 1), 'nnet.cnn.layer.ImageInputLayer')
            xRand = dlarray(xRand, "SSCB");
        elseif isa(net.Layers(1, 1), 'nnet.cnn.layer.FeatureInputLayer') || isa(net.Layers(1, 1), 'nnet.onnx.layer.FeatureInputLayer')
            xRand = dlarray(xRand, "CB");
        else
            disp(net.Layers(1,1));
            error("Unknown input format");
        end
    end
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

    precision = '%.16g'; % set the precision for all variables written to txt file
    % open file and start writing counterexamples
    fid = fopen(outputfile, 'a+');
    x = counterEx{1};
    x = reshape(x, [], 1);
    % begin specifying value for input example
    fprintf(fid,'(');
    for i = 1:length(x)
        fprintf(fid, "(X_" + string(i-1) + " " + num2str(x(i), precision)+ ")\n");
    end
    y = counterEx{2};
    y = reshape(y, [], 1);
    % specify values for output example
    for j =1:length(y)
        fprintf(fid, "(Y_" + string(j-1) + " " + num2str(y(j), precision)+ ")\n");
    end
    fprintf(fid, ')');
    % close and save file
    fclose(fid);

end

function counterEx = falsify_single(net, lb, ub, inputSize, nRand, Hs)
    counterEx = nan;
    xRand = create_random_examples(net, lb, ub, nRand, inputSize);
    s = size(xRand);
    n = length(s);
    %  look for counterexamples
    for i=1:s(n)
        x = get_example(xRand, i);
        yPred = predict(net, x);
        if isa(x, 'dlarray') % if net is a dlnetwork
            x = extractdata(x);
            yPred = extractdata(yPred);
        end
        % check if property violated
        yPred = reshape(yPred, [], 1); % convert to column vector (if needed)
        for h=1:length(Hs)
            if Hs(h).contains(double(yPred)) % property violated
                counterEx = {x; yPred}; % save input/output of countex-example
                break;
            end
        end
    end
    
end

function x = get_example(xRand,i)
    s = size(xRand);
    n = length(s);
    if n == 4
        x = xRand(:,:,:,i);
    elseif n == 3
        x = xRand(:,:,i);
    elseif n == 2
        x = xRand(:,i);
        xsize = size(x);
        if xsize(1) ~= 1 && ~isa(x,"dlarray")
            x = x';
        end
    else
        error("InputSize = "+string(s));
    end
end

