% Load instances
instances = readmatrix('instances.csv', 'OutputType','string');


for i=1:height(instances)

    status = -3;
    
    % Load network
    net = importNetworkFromONNX(instances(i,1), "InputDataFormats", "BC");
    nnvnet = matlab2nnv(net);
    
    % Load specification
    property = load_vnnlib(instances(i,2));
    lb = property.lb; % input lower bounds
    ub = property.ub; % input upper bounds
    prop = property.prop; % output spec to verify

    % Counterexample search
    t = tic;
    nR = 100; % number of samples
    IB = Box(lb,ub);
    xRand = IB.sample(nR-2);
    xRand = [lb, ub, xRand];
    xRand = dlarray(xRand, "CB");
    for k=1:nR
        x = xRand(:,k);
        yPred = predict(net, x);
        yPred = extractdata(yPred);
        yPred = reshape(yPred, [], 1); % convert to column vector (if needed)
        Hs = prop{1}.Hg;
        if Hs.contains(double(yPred)) % property violated
            status = 0;
            counterEx = {extractdata(x); yPred};
            break;
        end
    end
    
    t = toc(t);
    % disp("Counterexample search finished in "+string(t)+" seconds");

    % if status == 0
        % disp("-------     Counterexample found");
        % disp(counterEx{1});
        % disp(counterEx{2});
    % end
    
    % Create Input Set
    IS = Star(lb,ub);

    % Define reachability options
    reachOptions = struct;
    % reachOptions.reachMethod = 'relax-star-area';
    % reachOptions.relaxFactor = 0.0;
    % reachOptions.reachMethod = 'exact-star';
    reachOptions.reachMethod = 'approx-star';

    
    % Compute reachability
    t = tic;
    if status~=0
        ySet = nnvnet.reach(IS, reachOptions);
    end
    t = toc(t);
    
    % Verify property
    if status~=0
        status = verify_specification(ySet, prop);
    end

    if status == 2
        reachOptions = struct;
        reachOptions.reachMethod = 'exact-star';
        reachOptions.numCores = 4;
        ySet = nnvnet.reach(IS, reachOptions);
        status = verify_specification(ySet, prop);
    end


    % Print results
    disp("Instance ("+string(i)+") finished!! Result: "+string(status)+" , verification time: "+string(t));
    % disp("========================")

end

% Notes:
% some counterexamples
% most unknowns with relax-star (1)
% 