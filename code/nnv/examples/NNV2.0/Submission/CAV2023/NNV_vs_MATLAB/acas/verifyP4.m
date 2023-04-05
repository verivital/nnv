function verifyP4(reachOptionsList)

    % Setup the property
    resMAT = [];
    resNNV = [];
    timeMAT = zeros(45,1);
    n = length(reachOptionsList);
    timeNNV = zeros(45,n);
    
    XLower = [-0.303531156; -0.009549297; 0.0; 0.318181818; 0.083333333];
    XUpper = [-0.298552812; 0.009549297; 0.0; 0.5; 0.166666667];
    
    % NNV
    IS = ImageStar(XLower, XUpper);
    H = [1 -1 0 0 0; 1 0 -1 0 0; 1 0 0 -1 0; 1 0 0 0 -1];
    g = [0;0;0;0];
    
    % MATLAB
    XLower = dlarray(XLower, "CB");
    XUpper = dlarray(XUpper, "CB");
    
    % Iterate through all the networks to verify
    acasFolder = "onnx/";
    networks = dir(acasFolder);
    
    for i = 3:length(networks)
        % Load Network
        file = acasFolder + string(networks(i).name);
        net = importONNXNetwork(file, InputDataFormats='BCSS');

        % transform into NNV
        netNNV = matlab2nnv(net);

        % Transform the layers to fit MATLAB's algorithm
        Layers = net.Layers;
        Layers = Layers(5:end-1); % remove input and output layers
        input_layer = featureInputLayer(5); % number of inputs to acas xu
        elem_idxs = 2:3:length(Layers); % elementwiseLayers position
        for k=elem_idxs
            Layers(k-1).Bias = Layers(k).Offset; % add offset as bias on previous layer
        end
        Layers(elem_idxs) = []; % remove elementwiselayer
        netMAT = dlnetwork([input_layer; Layers]); % create dlnetwork for verification

        % Verify MATLAB
        t = tic;
        [YLower, YUpper] = estimateNetworkOutputBounds(netMAT, XLower, XUpper);
        res = verifyMAT(YLower, YUpper);
        timeMAT(i-2) = toc(t);
        resMAT = [resMAT; res];
        
        % Verify NNV
        res = [];
        for r = 1:n
            reachOpt = reachOptionsList(r);
            Rstar = [];
            t = tic;
            R = netNNV.reach(IS, reachOpt);
            for k = 1:length(R)
                Rstar = [Rstar R(k).toStar];
            end
            res = [res, verifyNNV(Rstar, H, g)];
            timeNNV(i-2,r) = toc(t);
        end
        resNNV = [resNNV; res];
        
    end

    % Save results
    save("results_p4", "resMAT", "resNNV", "timeMAT", "timeNNV", "reachOptionsList");

end

%% Helper Function

function result = verifyMAT(YLower, YUpper)
    if YLower(1) > YUpper(2) || YLower(1) > YUpper(3) || YLower(1) > YUpper(4) || YLower(1) > YUpper(5)
        result = categorical("violated"); % violated = safe
    elseif YUpper(1) <= YLower(2) && YUpper(1) <= YLower(3) && YUpper(1) <= YLower(4) && YUpper(1) <= YLower(5)
        result = categorical("verified"); % verified = unsafe
    else
        result = categorical("unproven"); % if approx methods used, then unproven, otherwise (exact) violated
    end
end

function result = verifyNNV(Set, H, b)
    for k = 1:length(Set)
        S = Set(k).intersectHalfSpace(H,b); % Violated?
        if isempty(S)
            result = categorical("verified"); % Does not interesct unsafe region (safe)
        elseif isempty(Set(k).intersectHalfSpace(-H,-b))
            result = categorical("violated");
            break;
        else
            result = categorical("nope");
            break;
        end
    end
end
