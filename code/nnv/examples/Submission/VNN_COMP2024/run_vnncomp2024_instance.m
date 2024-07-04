function [status, tTime] = run_vnncomp2024_instance(category, onnx, vnnlib, outputfile)

% single script to run all instances (supported) from the vnncomp2024

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

[net, nnvnet, needReshape] = load_vnncomp_network(category, onnx);

inputSize = net.Layers(1, 1).InputSize;


% Load property to verify
property = load_vnnlib(vnnlib);
lb = property.lb; % input lower bounds
ub = property.ub; % input upper bounds
prop = property.prop; % output spec to verify


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

% Define reachability options
% Let's try to choose this better from the get-go

% Option 1
% reachOptions_relax100 = struct;
% reachOptions_relax100.reachMethod = 'relax-star-range';
% reachOptions_relax100.relaxFactor = 1; 
% 
% % Option 2
% reachOptions_relax50 = struct;
% reachOptions_relax50.reachMethod = 'relax-star-range';
% reachOptions_relax50.relaxFactor = 0.5; 
% 
% % Option 3
% reachOptions_exact = struct;
% reachOptions_exact.reachMethod = 'exact-star';
% reachOptions_exact.reachOption = 'parallel';
% reachOptions_exact.numCores = feature('numcores');
% 
% % Option 4
% reachOptions_approx.reachMethod = 'approx-star';
% 
% % Choosing reachOptions (based on size, but not sure how to decice really...)
% if prod(inputSize) > 3000 % [32 32 3]
%     reachOptions = {reachOptions_relax100; reachOptions_relax50; reachOptions_approx};
% else
%     reachOptions = {reachOptions_relax50; reachOptions_approx; reachOptions_exact};
% end

reachOptions = struct;
reachOptions.lp_solver = "gurobi";
reachOptions.reachMethod = 'approx-star';
% reachOptions.reachMethod = 'exact-star';
% reachOptions.device = 'cpu';
% reachOptions.numCores = 24; % logical cores


% Check if property was violated earlier
if iscell(counterEx)
    status = 0;
end

vT = tic;

if status == 2 && isa(nnvnet, "NN") % no counterexample found and supported for reachability (otherwise, skip step 3 and write results)

% Choose how to verify based on vnnlib file
    if ~isa(lb, "cell") && length(prop) == 1 % one input, one output 

        IS = create_input_set(lb, ub, inputSize, needReshape);
        
        % Compute reachability
        ySet = nnvnet.reach(IS, reachOptions);

        % Verify property
        status = verify_specification(ySet, prop);

    elseif isa(lb, "cell") && length(lb) == length(prop) % multiple inputs, multiple outputs
        
        local_status = ones(length(lb),1); % track status for each specification in the vnnlib
        
        parfor spc = 1:length(lb) % We can compute these in parallel for faster computation

            lb_spc = lb{spc};
            ub_spc = ub{spc};
            
            IS = create_input_set(lb_spc, ub_spc, inputSize, needReshape);
    
            % Compute reachability
            ySet = nnvnet.reach(IS, reachOptions);
    
            % Verify property
            if isempty(ySet.C)
                dd = ySet.V; DD = ySet.V;
                ySet = Star(dd,DD);
            end

            % Add verification status
            local_status(spc) = verify_specification(ySet, prop(spc));

        end

        % Check for the global verification result
        if all(local_status == 1)
            status = 1;
        else
            status = 2;
        end

    elseif isa(lb, "cell") && length(prop) == 1 % one specification, multiple input definitions 

        local_status = ones(length(lb),1); % track status for each specification in the vnnlib
        
        parfor spc = 1:length(lb) % We can compute these in parallel for faster computation
            
            lb_spc = lb{spc};
            ub_spc = ub{spc};
            
            IS = create_input_set(lb_spc, ub_spc, inputSize, needReshape);
    
            % Compute reachability
            ySet = nnvnet.reach(IS, reachOptions);

            % Verify property
            local_status(spc) = verify_specification(ySet, prop);

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

if status == 2 && strcmp(reachOptions.reachMethod, 'exact-star')
    status = 0;
end

disp("Verification result: " + string(status));
disp("Counterexample search time: " + string(cEX_time));
disp("Reachability time: " + string(vT));
disp("Total Time: "+ string(tTime));
disp( " ");

% Write results to output file
if status == 0
    fid = fopen(outputfile, 'w');
    fprintf(fid, 'sat \n');
    fprintf(fid, '%f \n', tTime); % remove this line when running on submission site
    fclose(fid);
    write_counterexample(outputfile, counterEx)
elseif status == 1
    fid = fopen(outputfile, 'w');
    fprintf(fid, 'unsat \n');
    fprintf(fid, '%f \n', tTime); % remove this line when running on submission site
    fclose(fid);
elseif status == 2
    fid = fopen(outputfile, 'w');
    fprintf(fid, 'unknown \n');
    fprintf(fid, '%f \n', tTime); % remove this line when running on submission site
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
    IS = ImageStar(lb, ub); % 

end

function [net,nnvnet,needReshape] = load_vnncomp_network(category, onnx)
% load participating vnncomp 2024 benchmark NNs 
%
% Regular Track Benchmarks
% - 2024
%   - safeNLP
%   - cora
%   - linearizeNN
%   - cifar100
%   - tinyimagenet
% - 2023
%   - nn4sys
%   - dist-shift
%   - acasxu
%   - cgan
%   - collins_rul_cnn
%   - metaroom
%   - tllverifybench
%
% Extended Track Benchmarks
% - 2024
%   - ml4acopf
%   - dynaroars/benchmark-generation (name may change for competition)
%   - collins_yolo
%   - LSNC (not supported)
% - 2023
%   - yolo
%   - cctsdb_yolo
%   - collins_yolo_robustness (same as collins_yolo I believe)
%   - ml4acopf (different than this year's?)
%   - traffic_sign_recognition
%   - vggnet16
%   - vit   
%

    needReshape = 0; % default is to use MATLAB reshape, otherwise use the python reshape

    if contains(category, 'collins_rul')
        net = importNetworkFromONNX(onnx);
        nnvnet = matlab2nnv(net);
        if contains(onnx, 'full_window_40')
            needReshape = 2;
        else
            needReshape = 1;
        end

    elseif contains(category, "nn4sys")
        % nn4sys: onnx to matlab:
        if contains(onnx, "lindex")
            % nn4sys: onnx to nnv
            net = importNetworkFromONNX(onnx, "OutputDataFormats", "BC"); % lindex
            nnvnet = matlab2nnv(net);
        else
            error("We don't have those");
        end

    elseif contains(category, "ml4acopf")
        % ml4acopf: onnx to matlab
        net = importNetworkFromONNX(onnx, "InputDataFormats", "BC");
        nnvnet = "";

    elseif contains(category, "dist_shift")
        % dist_shift: onnx to matlab, , matlab to nnv?
        net = importNetworkFromONNX(onnx, "InputDataFormats", "BC", 'OutputDataFormats',"BC"); %reshape
        try 
            nnvnet = matlab2nnv(net);
        catch
            nnvnet = "";
        end

    elseif contains(category, "cgan")
        % cgan: onnx to nnv
        net = importNetworkFromONNX(onnx,"InputDataFormats", "BC", 'OutputDataFormats',"BC"); %reshape
        nnvnet = matlab2nnv(net);

    elseif contains(category, "vggnet16")
        % vgg16: onnx to matlab
        net = importNetworkFromONNX(onnx); % flattenlayer
        nnvnet = "";
        needReshape = 1;

    elseif contains(category, "tllverify")
        % tllverify: onnx to nnv
        net = importNetworkFromONNX(onnx,"InputDataFormats", "BC", 'OutputDataFormats',"BC");
        nnvnet = matlab2nnv(net);

    elseif contains(category, "vit")
        % vit: onnx to matlab
        net = importNetworkFromONNX(onnx, "InputDataFormats", "BCSS", 'OutputDataFormats',"BC");
        nnvnet = "";
        needReshape= 1;

    elseif contains(category, "cctsdb_yolo")
        % cctsdb_yolo: onnx to matlab
        error("We do not support this one");
        % net = importNetworkFromONNX(onnx);
        % nnvnet = "";
        % needReshape = ?;
        % We cannot support this one

    elseif contains(category, "collins_yolo")
        % collins_yolo: onnx to matlab
        net = importNetworkFromONNX(onnx);
        nnvnet = "";

    elseif contains(category, "yolo")
        % yolo: onnx to nnv
        net = importNetworkFromONNX(onnx); % padlayer
        nnvnet = matlab2nnv(net);

    elseif contains(category, "acasxu")
        % acasxu: onnx to nnv
        net = importNetworkFromONNX(onnx, "InputDataFormats","BCSS");
        nnvnet = matlab2nnv(net);

    elseif contains(category, "cifar100")
        % cifar100: onnx to nnv
        net = importNetworkFromONNX(onnx, "InputDataFormats","BCSS", "OutputDataFormats","BC");
        nnvnet = matlab2nnv(net);
        needReshape = 1;

    elseif contains(category, "tinyimagenet")
        % tinyimagenet: onnx to nnv
        net = importNetworkFromONNX(onnx, "InputDataFormats","BCSS", "OutputDataFormats","BC");
        nnvnet = matlab2nnv(net);
        % needReshpae = ?

    elseif contains(category, "LinearizeNN")
        % LinerizeNN: onnx to nnv
        net = importNetworkFromONNX(onnx, "InputDataFormats","BC", "OutputDataFormats","BC");
        nnvnet = matlab2nnv(net);
        % needReshape = ?

    elseif contains(category, "safeNLP")
        % safeNLP: onnx to nnv
        net = importNetworkFromONNX(onnx, "InputDataFormats","BC", "OutputDataFormats","BC");
        nnvnet = matlab2nnv(net);
        % needReshape = ?

    elseif contains(category, "cora")
        % cora benchmark: onnx 2 nnv
        net = importNetworkFromONNX(onnx, "InputDataFormats","BC", "OutputDataFormats","BC");
        nnvnet = matlab2nnv(net);
        % needReshape = ?

    elseif contains(category, "LSNC")
        % lyapunov benchmark: onnx to nnv (barely, some IR and opset version differences)
        net = importNetworkFromONNX(onnx, "InputDataFormats","BC", "OutputDataFormats","BC");
        nnvnet = "";
        % needReshape = ?

    elseif contains(category, "generation")
        % dynaroars/vnncomp-benchmark-generation
        net = importNetworkFromONNX(onnx, "InputDataFormats","BCSS", "OutputDataFormats","BC");
        nnvnet = matlab2nnv(net);
        % needReshape = ?

    elseif contains(category, "metaroom")
        % metaroom: onnx to matlab
        net = importNetworkFromONNX(onnx, "InputDataFormats","BCSS", "OutputDataFormats","BC");
        nnvnet = matlab2nnv(net);
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
        if needReshape ==2 % for collins only (full_window_40) and metaroom
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
        yPred = reshape(yPred, [], 1); % convert to column vector (if needed)
        % disp([x;yPred']);
        for h=1:length(Hs)
            if Hs(h).contains(double(yPred)) % property violated
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

