function [status, tTime] = run_vnncomp2025(category, onnx, vnnlib, outputfile)

% single script to run all instances (supported) from the vnncomp2024

t = tic; % start timer
status = 2; % unknown (to start with)

disp("We are running...")

% Process:
%  1) Load components
%  2) SAT? - Randomly evaluate network to search for counterexample 
%  3) UNSAT? 
%     a) Compute reachability 
%     b) Verify (compute intersection between reach set and output property halfspace)

%% 1) Load components

% Load networks

[net, nnvnet, needReshape, reachOptionsList, Prob] = load_vnncomp_network(category, onnx, vnnlib);

inputSize = net.Layers(1, 1).InputSize;

% Load property to verify
property = load_vnnlib(vnnlib);
lb = property.lb; % input lower bounds
ub = property.ub; % input upper bounds

prop = property.prop; % output spec to verify

                % reachOptions = reachOptionsList{1};
                % 
                % IS = create_input_set(lb, ub, inputSize, needReshape);
                % 
                % ySet = Prob_reach(nnvnet, IS, reachOptions);


%% 2) SAT?

nRand = 100; % number of random inputs (this can be changed)
% We got some penalties last year, why?
% Wrong vnnlib parsing? Wrong counterrexample writing? Do we need to reshape it?
% Let's test last years properties and make sure those errors/bugs are
% fixed before this year's submission

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
    warning("Working on adding support to other vnnlib properties");
end

cEX_time = toc(t);

%% 3) UNSAT?

% Check if property was violated earlier
if iscell(counterEx)
    status = 0;
end

vT = tic;

if status == 2  % no counterexample found and supported for reachability (otherwise, skip step 3 and write results)

% Choose how to verify based on vnnlib file
    if ~isa(lb, "cell") && length(prop) == 1 % one input, one output 

        if ~nnz(lb-ub) % lb == ub, not a set

            status = 1; % verified, since  we already tested this
            
        else

            while ~isempty(reachOptionsList)
                
                reachOptions = reachOptionsList{1};
    
                IS = create_input_set(lb, ub, inputSize, needReshape);
            
                % Compute reachability
                if Prob == 0
                    ySet = nnvnet.reach(IS, reachOptions);
                elseif Prob == 1
                    ySet = Prob_reach(nnvnet, IS, reachOptions);
                end
                % Verify property
                status = verify_specification(ySet, prop);
    
                if status == 1 % verified, then stop
                    break
                else
                    reachOptionsList = reachOptionsList(2:end);
                end
    
            end

        end

    elseif isa(lb, "cell") && length(lb) == length(prop) % multiple inputs, multiple outputs
        
        local_status = 2*ones(length(lb),1); % track status for each specification in the vnnlib
        
        parfor spc = 1:length(lb) % We can compute these in parallel for faster computation

            lb_spc = lb{spc};
            ub_spc = ub{spc};

            if ~nnz(lb_spc-ub_spc) % lb == ub, not a set

                local_status(spc) = 1; % verified, since we already tested this
                
            else

                reachOptPar = reachOptionsList;
                
                while ~isempty(reachOptPar)
    
                    reachOptions = reachOptPar{1};
                
                    IS = create_input_set(lb_spc, ub_spc, inputSize, needReshape);

                    % Compute reachability
                    if Prob == 0
                        ySet = nnvnet.reach(IS, reachOptions);
                    elseif Prob == 1
                        ySet = Prob_reach(nnvnet, IS, reachOptions);
                    end
            
                    % Verify property
                    if isempty(ySet.C)
                        dd = ySet.V; DD = ySet.V;
                        ySet = Star(dd,DD);
                    end
        
                    % Add verification status
                    tempStatus = verify_specification(ySet, prop(spc));
    
                    if tempStatus ~= 2 % verified, then stop (or falsified)
                        break
                    else
                        reachOptPar = reachOptPar(2:end);
                    end
    
                end

                local_status(spc) = tempStatus;

            end
            
        end

        % Check for the global verification result
        if all(local_status == 1)
            status = 1;
        else
            status = 2;
        end

    elseif isa(lb, "cell") && length(prop) == 1 % one specification, multiple input definitions 

        local_status = 2*ones(length(lb),1); % track status for each specification in the vnnlib, initialize as unknown
        
        parfor spc = 1:length(lb) % We can compute these in parallel for faster computation

            reachOptPar = reachOptionsList;
            
            lb_spc = lb{spc};
            ub_spc = ub{spc};

            if ~nnz(lb_spc-ub_spc) % lb == ub, not a set

                local_status(spc) = 1; % verified, since we already tested this earlier
                
            else
            
                while ~isempty(reachOptPar)
    
                    reachOptions = reachOptPar{1};
                
                    IS = create_input_set(lb_spc, ub_spc, inputSize, needReshape);

                    % Compute reachability
                    if Prob == 0
                        ySet = nnvnet.reach(IS, reachOptions);
                    elseif Prob == 1
                        ySet = Prob_reach(nnvnet, IS, reachOptions);
                    end
        
                    % Add verification status
                    tempStatus = verify_specification(ySet, prop(spc));
    
                    if tempStatus ~= 2 % verified, then stop (or falsified)
                        break
                    else
                        reachOptPar = reachOptPar(2:end);
                    end
    
                    local_status(spc) = tempStatus;
    
                end

            end

        end

        % Check for the global verification result
        if all(local_status == 1)
            status = 1;
        else
            status = 2;
        end

    else
        warning("Working on adding support to other vnnlib properties")
    end

end

vT = toc(vT);

%% 4) Process results

tTime = toc(t); % save total computation time

% if status == 2 && strcmp(reachOptions.reachMethod, 'exact-star')
%     status = 0;
% end

disp("Verification result: " + string(status));
disp("Counterexample search time: " + string(cEX_time));
disp("Reachability time: " + string(vT));
disp("Total Time: "+ string(tTime));
disp( " ");

% Write results to output file
if status == 0
    fid = fopen(outputfile, 'w');
    fprintf(fid, 'sat \n');
    % fprintf(fid, '%f \n', tTime); % remove this line when running on submission site
    fclose(fid);
    write_counterexample(outputfile, counterEx)
elseif status == 1
    fid = fopen(outputfile, 'w');
    fprintf(fid, 'unsat \n');
    % fprintf(fid, '%f \n', tTime); % remove this line when running on submission site
    fclose(fid);
elseif status == 2
    fid = fopen(outputfile, 'w');
    fprintf(fid, 'unknown \n');
    % fprintf(fid, '%f \n', tTime); % remove this line when running on submission site
    fclose(fid);
end

% quit; % does this work when running matlab.engine from python in background?


end

%% Helper functions

function IS = create_input_set(lb, ub, inputSize, needReshape)

    % Get input bounds
    if ~isscalar(inputSize)
        lb = reshape(lb, inputSize);
        ub = reshape(ub, inputSize);
    end

    % Format bounds into correct dimensions
    if needReshape == 1
        lb = permute(lb, [2 1 3]);
        ub = permute(ub, [2 1 3]);
    elseif needReshape == 2
        newSize = [inputSize(2) inputSize(1) inputSize(3:end)];
        lb = reshape(lb, newSize);
        lb = permute(lb, [2 1 3 4]);
        ub = reshape(ub, newSize);
        ub = permute(ub, [2 1 3 4]);
    end

    % Create input set
    
    equals = (lb==ub);

    if any(equals)
        diff = ub-lb;
        d = 10^(-6)*min(diff(:));
        lb(equals) = lb(equals) -d;
        ub(equals) = ub(equals) +d;
    end


    IS = ImageStar(lb, ub); 

end

function [net,nnvnet,needReshape,reachOptionsList, Prob] = load_vnncomp_network(category, onnx, vnnlib)



    needReshape = 0; % default is to use MATLAB reshape, otherwise use the python reshape

    numCores = feature('numcores'); % in case we select exact method

    if contains(category, "cgan")
        % cgan: onnx to nnv
        if ~contains(onnx, 'transformer')
            Prob = 0;
            net = importNetworkFromONNX(onnx,"InputDataFormats","BC","OutputDataFormats","BC");
            nnvnet = matlab2nnv(net);
            reachOptions.reachMethod = 'relax-star-area';
            reachOptions.relaxFactor = 0.8;
            reachOptionsList{1} = reachOptions;
            reachOptions = struct;
            reachOptions.reachMethod = 'approx-star'; % default parameters
            reachOptionsList{2} = reachOptions;            
        else
            Prob = 1;
            net = importNetworkFromONNX(onnx,"InputDataFormats","BC");
            nnvnet = net;
            reachOptions.train_epochs = 100;
            reachOptions.train_lr = 0.001;
            reachOptions.dims = [-1 -1];
            reachOptions.coverage = 0.999;
            reachOptions.confidence = 0.999;
            reachOptions.train_mode = 'Linear';
            reachOptions.surrogate_dim = [10, 10];
            reachOptions.threshold_normal = 1e-5;
            reachOptions.dlarrayType = 'CB';
            reachOptionsList{1} = reachOptions;

        end
        


    elseif contains(category, 'collins_aerospace_benchmark')
        needReshape = 2;
        net = importNetworkFromONNX(onnx, "InputDataFormats","BCSS");

        Prob = 0;
        try
            nnvnet = matlab2nnv(net);
        catch ME
            warning("matlab2nnv failed: %s", char(ME.message));
            nnvnet = net;  % Optional: clear or default
            Prob = 1;
        end
        
        if Prob == 0
            reachOptions = struct;
            reachOptions.reachMethod = 'approx-star'; % default parameters
            reachOptionsList{1} = reachOptions;
        elseif Prob == 1
            reachOptions.train_epochs = 100;
            reachOptions.train_lr = 0.001;
            reachOptions.coverage = 0.999;
            reachOptions.confidence = 0.999;
            reachOptions.train_mode = 'Linear';
            reachOptions.surrogate_dim = [];
            reachOptions.threshold_normal = 1e-5;

            if isa(net.Layers(1, 1), 'nnet.cnn.layer.ImageInputLayer')
                reachOptions.dlarrayType = 'SSCB';
            elseif isa(net.Layers(1, 1), 'nnet.cnn.layer.FeatureInputLayer') || isa(net.Layers(1, 1), 'nnet.onnx.layer.FeatureInputLayer')
                reachOptions.dlarrayType = 'CB';
            else
                disp(net.Layers(1,1));
                error("Unknown input format");
            end

            reachOptionsList{1} = reachOptions;
        end


    elseif contains(category, "LCNS_ReLU")
        % LCNS_ReLU: onnx to matlab:

        net = importNetworkFromONNX(onnx, "InputDataFormats","BC", "OutputDataFormats","BC");

        Prob = 0;
        try
            nnvnet = matlab2nnv(net);
        catch ME
            warning("matlab2nnv failed: %s", char(ME.message));
            nnvnet = net;  % Optional: clear or default
            Prob = 1;
        end
        
        if Prob == 0
            reachOptions = struct;
            reachOptions.reachMethod = 'approx-star';
            reachOptionsList{1} = reachOptions;
        elseif Prob == 1
            reachOptions.train_epochs = 100;
            reachOptions.train_lr = 0.001;
            reachOptions.coverage = 0.999;
            reachOptions.confidence = 0.999;
            reachOptions.train_mode = 'Linear';
            reachOptions.surrogate_dim = [];
            reachOptions.threshold_normal = 1e-5;

            if isa(net.Layers(1, 1), 'nnet.cnn.layer.ImageInputLayer')
                reachOptions.dlarrayType = 'SSCB';
            elseif isa(net.Layers(1, 1), 'nnet.cnn.layer.FeatureInputLayer') || isa(net.Layers(1, 1), 'nnet.onnx.layer.FeatureInputLayer')
                reachOptions.dlarrayType = 'CB';
            else
                disp(net.Layers(1,1));
                error("Unknown input format");
            end

            reachOptionsList{1} = reachOptions;
        end
        

    elseif contains(category, "ml4acopf")
        % ml4acopf: onnx to matlab

        net = importNetworkFromONNX(onnx, "InputDataFormats","BC", "OutputDataFormats","BC");

        Prob = 0;
        try
            nnvnet = matlab2nnv(net);
        catch ME
            warning("matlab2nnv failed: %s", char(ME.message));
            nnvnet = net;  % Optional: clear or default
            Prob = 1;
        end
        
        if Prob == 0
            reachOptions = struct;
            reachOptions.reachMethod = 'approx-star';
            reachOptionsList{1} = reachOptions;
        elseif Prob == 1
            reachOptions.train_epochs = 500;
            reachOptions.train_lr = 0.0001;
            reachOptions.coverage = 0.999;
            reachOptions.confidence = 0.999;
            reachOptions.train_mode = 'Linear';
            reachOptions.surrogate_dim = [-1,-1];
            reachOptions.threshold_normal = 1e-5;

            if isa(net.Layers(1, 1), 'nnet.cnn.layer.ImageInputLayer')
                reachOptions.dlarrayType = 'SSCB';
            elseif isa(net.Layers(1, 1), 'nnet.cnn.layer.FeatureInputLayer') || isa(net.Layers(1, 1), 'nnet.onnx.layer.FeatureInputLayer')
                reachOptions.dlarrayType = 'CB';
            else
                disp(net.Layers(1,1));
                error("Unknown input format");
            end

            reachOptionsList{1} = reachOptions;
        end


    elseif contains(category, "Sound")

        net = importNetworkFromONNX(onnx, "InputDataFormats","BC", "OutputDataFormats","BC");
        nnvnet = net;  % Optional: clear or default
        Prob = 1;
        reachOptions.train_epochs = 500;
        reachOptions.train_lr = 0.001;
        reachOptions.coverage = 0.999;
        reachOptions.confidence = 0.999;
        reachOptions.train_mode = 'Linear';
        reachOptions.surrogate_dim = [-1,-1];
        reachOptions.threshold_normal = 1e-5;

        if isa(net.Layers(1, 1), 'nnet.cnn.layer.ImageInputLayer')
            reachOptions.dlarrayType = 'SSCB';
        elseif isa(net.Layers(1, 1), 'nnet.cnn.layer.FeatureInputLayer') || isa(net.Layers(1, 1), 'nnet.onnx.layer.FeatureInputLayer')
            reachOptions.dlarrayType = 'CB';
        else
            disp(net.Layers(1,1));
            error("Unknown input format");
        end

        reachOptionsList{1} = reachOptions;

        % net = importNetworkFromONNX(onnx, "InputDataFormats","BC", "OutputDataFormats","BC");
        % nnvnet = matlab2nnv2(net);
        % reachOptions = struct;
        % reachOptions.reachMethod = 'approx-star';
        % reachOptionsList{1} = reachOptions;

    elseif contains(category, "tinyimagenet")
        % tinyimagenet: onnx to nnv
        net = importNetworkFromONNX(onnx, "InputDataFormats","BCSS", "OutputDataFormats","BC");
        nnvnet = net;  % Optional: clear or default
        Prob = 1;
        reachOptions.train_epochs = 150;
        reachOptions.train_lr = 0.0001;
        reachOptions.coverage = 0.999;
        reachOptions.confidence = 0.999;
        reachOptions.train_mode = 'Linear';
        reachOptions.surrogate_dim = [-1,-1];
        reachOptions.threshold_normal = 1e-5;

        if isa(net.Layers(1, 1), 'nnet.cnn.layer.ImageInputLayer')
            reachOptions.dlarrayType = 'SSCB';
        elseif isa(net.Layers(1, 1), 'nnet.cnn.layer.FeatureInputLayer') || isa(net.Layers(1, 1), 'nnet.onnx.layer.FeatureInputLayer')
            reachOptions.dlarrayType = 'CB';
        else
            disp(net.Layers(1,1));
            error("Unknown input format");
        end

        reachOptionsList{1} = reachOptions;
        needReshape = 2;
    
    elseif contains(category, "ViT")
        % vit: onnx to matlab
        needReshape = 2;
        net = importNetworkFromONNX(onnx, "InputDataFormats", "BCSS", 'OutputDataFormats',"BC");

        Prob = 0;
        try
            nnvnet = matlab2nnv(net);
        catch ME
            warning("matlab2nnv failed: %s", char(ME.message));
            nnvnet = net;  % Optional: clear or default
            Prob = 1;
        end
        
        if Prob == 0
            reachOptions = struct;
            reachOptions.reachMethod = 'approx-star';
            reachOptionsList{1} = reachOptions;
        elseif Prob == 1
            reachOptions.train_epochs = 100;
            reachOptions.train_lr = 0.001;
            reachOptions.coverage = 0.999;
            reachOptions.confidence = 0.999;
            reachOptions.train_mode = 'Linear';
            reachOptions.surrogate_dim = [-1,-1];
            reachOptions.threshold_normal = 1e-5;

            if isa(net.Layers(1, 1), 'nnet.cnn.layer.ImageInputLayer')
                reachOptions.dlarrayType = 'SSCB';
            elseif isa(net.Layers(1, 1), 'nnet.cnn.layer.FeatureInputLayer') || isa(net.Layers(1, 1), 'nnet.onnx.layer.FeatureInputLayer')
                reachOptions.dlarrayType = 'CB';
            else
                disp(net.Layers(1,1));
                error("Unknown input format");
            end

            reachOptionsList{1} = reachOptions;
        end

    elseif contains(category, "yolo")
        % yolo: onnx to nnv
        net = importNetworkFromONNX(onnx, "InputDataFormats","BCSS", "OutputDataFormats","BC");
        nnvnet = net;  % Optional: clear or default
        Prob = 1;
        reachOptions.train_epochs = 200;
        reachOptions.train_lr = 0.0001;
        reachOptions.coverage = 0.999;
        reachOptions.confidence = 0.999;
        reachOptions.train_mode = 'Linear';
        reachOptions.surrogate_dim = [-1,-1];
        reachOptions.threshold_normal = 1e-5;

        if isa(net.Layers(1, 1), 'nnet.cnn.layer.ImageInputLayer')
            reachOptions.dlarrayType = 'SSCB';
        elseif isa(net.Layers(1, 1), 'nnet.cnn.layer.FeatureInputLayer') || isa(net.Layers(1, 1), 'nnet.onnx.layer.FeatureInputLayer')
            reachOptions.dlarrayType = 'CB';
        else
            disp(net.Layers(1,1));
            error("Unknown input format");
        end

        reachOptionsList{1} = reachOptions;
        needReshape = 2;

    else % all other benchmarks
        % traffic: onnx to matlab: opset15 issues
        error("ONNX model not supported")
    end

end

% Create an array of random examples from input set and reshape if necessary
% We use dlnetwork for simulation (MATLAB data structure)
function xRand = create_random_examples(net, lb, ub, nR, inputSize, needReshape)
    xB = Box(lb, ub); % lb, ub must be vectors
    xRand = xB.sample(nR-2);
    xRand = [lb, ub, xRand];
    if needReshape
        if needReshape ==2 % for collins only (full_window_40) and metaroom  flattened version of CHW (in python)--> HWC in matlab
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

% Write counterexample to output file
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

% Falsification function (random simulation looking for counterexamples)
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
        yPred = reshape(yPred, [], 1); % convert to column vector (if needed) The specification is also matlab translated, so if you use toStar() then your reachset is applicable
        % disp([x;yPred']);
        for h=1:length(Hs)
            if Hs(h).contains(double(yPred)) % property violated
                % check if the counter example needs to be reshaped
                n = numel(x);
                if needReshape == 2
                    % x = reshape(x, [n 1]);
                    x = permute(x, [2 1 3]);
                elseif needReshape == 1
                    if ndims(x) == 3 % RGB  image
                        x = permute(x, [2 1 3]);
                    end
                end
                counterEx = {x; yPred}; % save input/output of countex-example
                break;
            end
        end
    end
    
end

% Get random example from input set
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

