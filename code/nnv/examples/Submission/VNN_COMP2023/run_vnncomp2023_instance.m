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

[net, nnvnet, needReshape] = load_vnncomp_network(category, onnx);
inputSize = net.Layers(1, 1).InputSize;

cd(old_path); % go back to where we ran the functions from

% Load property to verify
property = load_vnnlib(vnnlib);
lb = property.lb; % input lower bounds
ub = property.ub; % input upper bounds
prop = property.prop; % output spec to verify


%% 2) SAT?

nRand = 100; % number of random inputs

% Choose how to falsify based on vnnlib file
if ~isa(lb, "cell") && length(prop) == 1 % one input, one output 
    counterEx = falsify_single(net, lb, ub, inputSize, nRand, prop{1}.Hg, needReshape);
elseif isa(lb, "cell") && length(lb) == length(prop) % multiple inputs, multiple outputs
    for spc = 1:length(lb) % try parfeval, parfor does not work for early return
        counterEx = falsify_single(net, lb{spc}, ub{spc}, inputSize, nRand, prop{spc}.Hg, needReshape);
        if iscell(counterEx)
            break
        end
    end
elseif isa(lb, "cell") && length(prop) == 1 % can violate the output property from any of the input sets
    for arr = 1:length(lb) % try parfeval, parfor does not work for early return
        counterEx = falsify_single(net, lb{arr}, ub{arr}, inputSize, nRand, prop{1}.Hg, needReshape);
        if iscell(counterEx)
            break
        end
    end
else
    warning("Working on adding support to other vnnlib properties")
end

toc(t);

%% 3) UNSAT?

% Define reachability options
% Could have chosen a bit better based on benchmarks, but it is what it is

% Option 1
reachOptions_relax100 = struct;
reachOptions_relax100.reachMethod = 'relax-star-range';
reachOptions_relax100.relaxFactor = 1; 

% Option 2
reachOptions_relax50 = struct;
reachOptions_relax50.reachMethod = 'relax-star-range';
reachOptions_relax50.relaxFactor = 0.5; 

% Option 3
reachOptions_exact = struct;
reachOptions_exact.reachMethod = 'exact-star';
reachOptions_exact.reachOption = 'parallel';
reachOptions_exact.numCores = feature('numcores');

% Option 4
reachOptions_approx.reachMethod = 'approx-star';

% Choosing reachOptions (based on size, but not sure how to decice really...)
if prod(inputSize) > 3000 % [32 32 3]
    reachOptions = {reachOptions_relax100; reachOptions_relax50; reachOptions_approx};
else
    reachOptions = {reachOptions_relax50; reachOptions_approx; reachOptions_exact};
end


% Check if property was violated earlier
if iscell(counterEx)
    status = 0;
end

if status == 2 && isa(nnvnet, "NN") % no counterexample found and supported for reachability (otherwise, skip step 3 and write results)

% Choose how to verify based on vnnlib file
    if ~isa(lb, "cell") && length(prop) == 1 % one input, one output 
        % Get input set
        if ~isscalar(inputSize)
            lb = reshape(lb, inputSize);
            ub = reshape(ub, inputSize);
        end
        if needReshape
            lb = permute(lb, [2 1 3]);
            ub = permute(ub, [2 1 3]);
        end
        IS = ImageStar(lb, ub);
        % Compute reachability
        ySet = nnvnet.reach(IS, reachOptions{1});
        % Verify property
        status = verify_specification(ySet, prop);
        if status == 2 % refine if unknown
            % Compute reachability
            ySet = nnvnet.reach(IS, reachOptions{2});
            % Verify property
            status = verify_specification(ySet, prop);
            if status == 2 % refine if unknown
                % Compute reachability
                ySet = nnvnet.reach(IS, reachOptions{3});
                % Verify property
                status = verify_specification(ySet, prop);
            end
        end

    elseif isa(lb, "cell") && length(lb) == length(prop) % multiple inputs, multiple outputs
        local_status = ones(length(lb),1);
        parfor spc = 1:length(lb)
            % Get input set
            if ~isscalar(inputSize)
                lbS = reshape(lb{spc}, inputSize);
                ubS = reshape(ub{spc}, inputSize);
            else
                lbS = lb{spc};
                ubS = ub{spc};
            end
            if needReshape
                lbS = permute(lbS, [2 1 3]);
                ubS = permute(ubS, [2 1 3]);
            end
            IS = ImageStar(lbS, ubS);
            % Compute reachability
            ySet = nnvnet.reach(IS, reachOptions{1});
            % Verify property
            if isempty(ySet.C)
                dd = ySet.V; DD = ySet.V;
                ySet = Star(dd,DD);
            end
            local_status(spc) = verify_specification(ySet, prop(spc));
            if local_status(spc) == 2 % refine if unknown
                % Compute reachability
                ySet = nnvnet.reach(IS, reachOptions{2});
                % Verify property
                local_status(spc) = verify_specification(ySet, prop(spc));
                if local_status(spc) == 2 % refine if unknown
                    % Compute reachability
                    ySet = nnvnet.reach(IS, reachOptions{3});
                    % Verify property
                    local_status(spc) = verify_specification(ySet, prop(spc));
                end
            end
        end
        if all(local_status == 1)
            status = 1;
        else
            status = 2;
        end
    elseif isa(lb, "cell") && length(prop) == 1
        local_status = ones(length(lb),1);
        parfor spc = 1:length(lb)
            % Get input set
            if ~isscalar(inputSize)
                lbS = reshape(lb{spc}, inputSize);
                ubS = reshape(ub{spc}, inputSize);
            else
                lbS = lb{spc};
                ubS = ub{spc};
            end
            if needReshape
                lbS = permute(lbS, [2 1 3]);
                ubS = permute(ubS, [2 1 3]);
            end
            IS = ImageStar(lbS, ubS);
            % Compute reachability
            ySet = nnvnet.reach(IS, reachOptions{1});
            % Verify property
            local_status(spc) = verify_specification(ySet, prop);
            if local_status(spc) == 2 % refine if unknown
                % Compute reachability
                ySet = nnvnet.reach(IS, reachOptions{2});
                % Verify property
                local_status(spc) = verify_specification(ySet, prop);
                if local_status(spc) == 2 % refine if unknown
                    % Compute reachability
                    ySet = nnvnet.reach(IS, reachOptions{3});
                    % Verify property
                    local_status(spc) = verify_specification(ySet, prop);
                end
            end
        end
        if all(local_status == 1)
            status = 1;
        else
            status = 2;
        end
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

% quit; % does this work when running matlab.engine from python in background?


end

function [net,nnvnet, needReshape] = load_vnncomp_network(category, onnxFile)
% load vnncomp 2023 and 2022 benchmark NNs (subset support)

    needReshape = 0; % default is to use MATLAB reshape, otherwise use the python reshape

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
        if contains(onnx, 'full_window_40')
            needReshape = 2;
        else
            needReshape = 1;
        end

    elseif contains(category, "ml4acopf")
        net = load(onnx);
        net = net.net;
        nnvnet = "";

    elseif contains(category, "nn4sys")
        % nn4sys: onnx to matlab
        if contains(onnxFile, "lindex")
            net = load(onnx);
            net = assembleNetwork(net.net);
            nnvnet = matlab2nnv(net);
        else
            error("We don't have those");
        end
        
    elseif contains(category, "dist_shift")
        % dist_shift: onnx to matlab:
        net = load(onnx);
        net = net.net;
        try 
            nnvnet = matlab2nnv(net);
        catch
            nnvnet = "";
        end
        
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
        try
            nnvnet = matlab2nnv(net);
        catch
            nnvnet = "";
        end
        needReshape = 1;
        
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
        needReshape= 1;
        
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

    elseif contains(category, 'sri_resnet')
        net = load(onnx);
        net = net.net;
        nnvnet = matlab2nnv(net);
        needReshape = 1;
        
    elseif contains(category, 'biasfield')
        % onnx to nnv
        if contains(onnx, "base")
            needReshape = 1;
        end
        net = load(onnx);
        net = net.net;
        nnvnet = matlab2nnv(net);

    elseif contains(category, 'rl_bench')
        % rl: onnx to nnv
        net = load(onnx);
        net = net.net;
        nnvnet = matlab2nnv(net);

    elseif contains(category, 'reach_prob')
        % reach_prob: onnx to nnv
        net = load(onnx);
        net = net.net;
        nnvnet = matlab2nnv(net);

    elseif contains(category, 'mnist_fc')
        % mnist: onnx to nnv (input is different here, need to check our
        % code still works for it)
        net = load(onnx);
        net = net.net;
        nnvnet = matlab2nnv(net);

    elseif contains(category, 'cifar100')
        net = load(onnx);
        net = net.net;
        nnvnet = matlab2nnv(net);
        needReshape = 1;

    elseif contains(category, 'oval')
        net = load(onnx);
        net = net.net;
        nnvnet = matlab2nnv(net);
        needReshape = 1;

    elseif contains(category, 'cifar2020')
        % onnx to nnv
        net = load(onnx);
        net = net.net;
        nnvnet = matlab2nnv(net);
        needReshape = 1;

    else % all other benchmarks
        % traffic: onnx to matlab: opset15 issues
        error("ONNX model not supported")
    end

end

% Had to change this because importer fails using the command line...s
% function [net,nnvnet] = load_vnncomp_network_local(category, onnx)
% % load vnncomp 2023 benchmark NNs (subset support)
% 
%     % collins_rul: onnx to nnvnet
%     % collins_nets = load_collins_NNs;
%     if contains(category, 'collins_rul')
%         net = importONNXNetwork(onnx);
%         nnvnet = matlab2nnv(net);
% 
%     elseif contains(category, "nn4sys")
%         % nn4sys: onnx to matlab:
%         net = importONNXLayers(onnx, "OutputDataFormats", "BC"); % lindex
%         nnvnet = matlab2nnv(net);
%         
%     elseif contains(category, "dist_shift")
%         % dist_shift: onnx to matlab:
%         net = importONNXNetwork(onnx, "InputDataFormats", "BC", 'OutputDataFormats',"BC"); %reshape
%         nnvnet = "";
%         
%     elseif contains(category, "cgan")
%         % cgan
%         net = importONNXNetwork(onnx,"InputDataFormats", "BC", 'OutputDataFormats',"BC"); %reshape
%         nnvnet = matlab2nnv(net);
%         
%     elseif contains(category, "vgg16")
%         % vgg16: onnx to matlab
%         net = importONNXNetwork('vgg16-7.onnx'); % flattenlayer
%         %reshapedInput = python_reshape(input,net_vgg.Layers(1,1).InputSize); % what is the input? assume it's all the same?
%         %nnvnet = matlab2nnv(net);
%         nnvnet = "";
%         
%     elseif contains(category, "tllverify")
%         % tllverify: onnx to matlab
%         net = importONNXNetwork(onnx,"InputDataFormats", "BC", 'OutputDataFormats',"BC");
%         nnvnet = matlab2nnv(net);
%         
%     elseif contains(category, "vit")
%         % vit: onnx to matlab
%         net = importONNXNetwork(onnx, "TargetNetwork","dlnetwork" );
%         nnvnet = "";
%         
%     elseif contains(category, "cctsdb_yolo")
%         % cctsdb_yolo: onnx to matnet
%         net = importONNXNetwork(onnx, "TargetNetwork","dlnetwork" );
%         nnvnet = "";
%         
%     elseif contains(category, "collins_yolo")
%         % collins_yolo: onnx to matlab:
%         net = importONNXNetwork(onnx, "TargetNetwork","dlnetwork" );
%         nnvnet = "";
% 
%     elseif contains(category, "yolo")
%         % yolo: onnx to matlab
%         net = importONNXNetwork(onnx); % padlayer
%         nnvnet = matlab2nnv(net);
% 
%     elseif contains(category, "acasxu")
%         % acasxu: onnx to nnv
%         net = importONNXNetwork(onnx, "InputDataFormats","BCSS");
%         nnvnet = matlab2nnv(net);
%         
%     else % all other benchmarks
%         % traffic: onnx to matlab: opset15 issues
%         error("ONNX model not supported")
%     end
% 
% end

function xRand = create_random_examples(net, lb, ub, nR, inputSize, needReshape)
    xB = Box(lb, ub); % lb, ub must be vectors
    xRand = xB.sample(nR-2);
    xRand = [lb, ub, xRand];
    if needReshape
        if needReshape ==2 % for collins only (full_window_40)
            newSize = [inputSize(2) inputSize(1) inputSize(3:end)];
            xRand = reshape(xRand, [newSize nR]);
            xRand = permute(xRand, [2 1 3 4]);
        else
            xRand = reshape(xRand, [inputSize nR]);
            xRand = permute(xRand, [2 1 3 4]);
        end
%         xRand = python_reshape(xRand, [inputSize nR]);
    else
        xRand = reshape(xRand,[inputSize nR]); % reshape vectors into net input size
    end
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

function counterEx = falsify_single(net, lb, ub, inputSize, nRand, Hs, needReshape)
    counterEx = nan;
    xRand = create_random_examples(net, lb, ub, nRand, inputSize, needReshape);
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

