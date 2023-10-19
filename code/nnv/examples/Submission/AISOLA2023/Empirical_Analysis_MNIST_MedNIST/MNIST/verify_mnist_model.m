function verify_mnist_model(modelpath, Sets, targets)

    % Get names for saving and display
    modelName = split(modelpath, filesep);
    regName = modelName{end-3};
    initName = modelName{end-2};
    saveName = split(modelName{end},'.');
    saveName = saveName{1};
    disp(['Verifying model with regularization: ', regName, ' , initialization: ' , initName, ', name: ', saveName]);
    
    % Load NN
    load(modelpath); % loads net + accuracy
    nn = matlab2nnv(net); % transform net to nnv format (NN)
    % ensure I/O are correct
    nn.InputSize = net.Layers(1).InputSize;
    if isa(net, "SeriesNetwork")
        nn.OutputSize = net.Layers(end-2).OutputSize; % dropout and l2
    else
        nn.OutputSize = net.Layers(end-1).OutputSize; % jacobian
    end

    % Check input data
    if length(Sets) ~= length(targets)
        error("Input and outputs must match");
    end
    
    % Initialize save vars
    N = length(targets);
    res = zeros(N,2); % col 1: result, col 2: time
    
    % Define reachability options
    reachOptions = struct;
    reachOptions.reachMethod = 'approx-star';
    
    % Begin verification
    for i=1:N
        if ~mod(i, 20)
            disp("Verifying input "+string(i)+ " ...");
        end
        try
            t = tic;
            res(i,1) = nn.verify_robustness(Sets(i), reachOptions, targets(i));
            res(i,2) = toc(t);
        catch
            res(i,1) = -1;
            res(i,2) = -1;
        end
    end
    
    % Save results
    save(['results', filesep, 'rob_', saveName], 'res');

    % Show results in command window
    disp("========  RESULTS  ========");
    disp("Model: "+string(saveName));
    disp("");
    disp("Average computation time: "+string(sum(res(:,2))/N));
    disp("Robust = "     + string(sum(res(:,1)==1)) + " out of " + string(N) + " images");
    disp("Unknown = "    + string(sum(res(:,1)==2)) + " out of " + string(N) + " images");
    disp("Not Robust = " + string(sum(res(:,1)==0)) + " out of " + string(N) + " images");
    disp(" ");

end
