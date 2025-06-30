% Load instances
instances = readmatrix('instances.csv', 'OutputType','string', "Delimiter",",");


for i=1:height(instances)

    status = -3;
    
    % Load network
    net = importNetworkFromONNX(instances(i,1), "InputDataFormats", "BCSS");
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
    xRand = dlarray(xRand, "SSCB");
    for k=1:nR
        x = xRand(:,k);
        x = reshape(x, [64,64,1,1]);
        x = permute(x, [2 1 3]);
        yPred = predict(net, x);
        yPred = extractdata(yPred);
        yPred = reshape(yPred, [], 1); % convert to column vector (if needed)
        Hs = prop{1}.Hg;
        for h = 1:length(Hs)
            if Hs(h).contains(double(yPred)) % property violated
                status = 0;
                counterEx = {extractdata(x); yPred};
                break;
            end
        end
    end
    
    t = toc(t);
    disp("Counterexample search finished in "+string(t)+" seconds");

    if status == 0
        "Counterexample found"
        % disp(counterEx{1});
        % disp(counterEx{2});
    end
    
    % Create Input Set
    lb = reshape(lb, [64,64]);
    lb = permute(lb, [2 1 3]);
    ub = reshape(ub, [64,64]);
    ub = permute(ub, [2 1 3]);
    IS = ImageStar(lb,ub);

    % Define reachability options
    reachOptions = struct;
    reachOptions.reachMethod = 'relax-star-area';
    reachOptions.relaxFactor = 1;
    % reachOptions.reachMethod = 'exact-star';
    
    % Compute reachability
    t = tic;
    ySet = nnvnet.reach(IS, reachOptions);
    t = toc(t);
    
    % Verify property
    if status~=0
        status = verify_specification(ySet, prop);
    end


    % Print results
    disp("Instance ("+string(i)+") finished!! Result: "+string(status)+" , verification time: "+string(t));
    disp("========================")

end

% Notes:
% Some counterexamples found
% Majority can be verified with relax-star-area (1)
% Can refine to approx-star after initial run
