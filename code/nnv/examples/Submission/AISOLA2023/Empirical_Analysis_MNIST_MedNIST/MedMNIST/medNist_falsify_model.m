function medNist_falsify_model(modelpath, xRand, targets)

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
    if length(xRand) ~= length(targets)
        error("Input and outputs must match");
    end
    
    % Initialize save vars
    N = length(targets);
    res = zeros(N,2); % col 1: result, col 2: time
    counterExamples = cell(N,1);
    
    % Begin verification
    for i=1:N
        if ~mod(i, 20)
            disp("Verifying input "+string(i)+ " ...");
        end
        t = tic;
        cEx = falsify_single(nn, xRand{i}, targets(i));
        if isnan(cEx{1})
            res(i,1) = 2; % counterExample not found
        else
            res(i,1) = 0; % counterExample found
        end
        res(i,2) = toc(t);
        counterExamples(i) = cEx;
    end
    
    % Save results
    save(['results_falsify', filesep, 'falsify_', saveName], 'res', 'counterExamples');

end

function counterEx = falsify_single(net, xImg, yTarget)
    counterEx = {nan};
    s = size(xImg);
    %  look for counterexamples
    for i=1:s(2)
        x = xImg(:,1);
        x = reshape(x, [64, 64]);
        yPred = net.evaluate(x);
        % check if property violated
        [~,yPred] = max(yPred);
        if yPred ~= yTarget
            counterEx = {x}; % save input/output of countex-example
            break;
        end
    end
end
