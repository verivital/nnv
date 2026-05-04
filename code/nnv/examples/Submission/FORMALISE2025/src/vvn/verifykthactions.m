function [res, time, met] = verifykthactions(smpLen, verAlg, index, epsIndex)

    if smpLen ~= 8 && smpLen ~= 16 && smpLen ~= 32
        printf("smpLen argument was invalid. Must be 8, 16, or 32.")
        return
    end
    
    if verAlg ~= "relax" && verAlg ~= "approx"
        printf("verAlg argument was invalid. Must be 'relax' or 'approx'.")
        return
    end

    fprintf("Running robustness verification on KTH Actions %df dataset...", smpLen);

    % Load data
    data = readNPY(sprintf("data/KTHActions/kthactions_grayscale_%df_verification_data.npy", smpLen));
    labels = readNPY(sprintf("data/KTHActions/kthactions_grayscale_%df_verification_labels.npy", smpLen));

    % Experimental variables
    numClasses = 6;

    % Size of attack
    epsilon = [1/255; 2/255; 3/255];
    nE = length(epsilon);

    % Load the model
    modelName = sprintf("kthactions_c3d_%df.onnx", smpLen);
    netonnx = importONNXNetwork("models/" + modelName, "InputDataFormats", "TBCSS", "OutputDataFormats", "BC");
    net = matlab2nnv(netonnx);
    net.OutputSize = numClasses;
    disp("Finished loading model: " + modelName);

    % Get the sample
    s = data(index,:,:,:,:);
    s = squeeze(s);
    label = labels(index) + 1;

    output = net.evaluate(s);
    [~, P] = max(output);

    % choose the epsilon vaue
    eps = epsilon(epsIndex);

    % Perform L_inf attack
    VS = L_inf_attack(s, eps, smpLen);

    %%%%%%%%%%%%%%%%%
    % FALSIFICATION %
    %%%%%%%%%%%%%%%%%

    t = tic;

    % initialize return variables
    met = verAlg;
    res = 100;
    time = 0.0;

    % convert the model to dlnetwork and the data to dlarray
    lgraph = layerGraph(netonnx);
    lastL = lgraph.Layers(end);
    lgraph = removeLayers(lgraph, lastL.Name);
    dlnet = dlnetwork(lgraph);

    s_pgd = reshape(s, 1, 1, smpLen, 120, 160); % prepare sample for PGD
    s_pgd = dlarray(s_pgd, 'BCSSS');

    T = dlarray(onehotencode(label, 1, 'ClassNames', 1:numClasses)); % one-hot encoded label vector

    % generate adversarial sample
    Xadv = PGD(dlnet, s_pgd, T, ...
        'Epsilon', eps, 'Alpha', eps / 4, 'Steps', 100, ...
        'DataMin', 0, 'DataMax', 1);
    
    advOutput = net.evaluate(Xadv);
    [~, advP] = max(advOutput);

    if advP ~= label
        res = 0;
        time = toc(t);
        met = 'falsified';
        return;
    end

    %%%%%%%%%%%%%%%%
    % VERIFICATION %
    %%%%%%%%%%%%%%%%

    % Verification settings
    reachOptions = struct;
    if verAlg == "relax"
        reachOptions.reachMethod = "relax-star-area";
        reachOptions.relaxFactor = 0.5;
    elseif verAlg == "approx"
        reachOptions.reachMethod = "approx-star";
    end

    fprintf('Starting reachability analysis with epsilon %d \n', eps);
    
    try
        % run verification algorithm
        temp = net.verify_robustness(VS, reachOptions, labels(index)+1);
                
    catch ME
        met = ME.message;
        temp = -1;
    end
    
    res = temp;
    time = toc(t);

end

%% Helper Functions
function VS = L_inf_attack(x, epsilon, numFrames)
    lb = squeeze(x);
    ub = squeeze(x);

    % Perturb the frames
    for fn=1:numFrames
        lb(fn, :, :, :) = x(fn, :, :, :) - epsilon;
        ub(fn, :, :, :) = x(fn, :, :, :) + epsilon;
    end

    % Clip the perturbed values to be between 0-1
    lb_min = zeros(numFrames, 120, 160, 1);
    ub_max = ones(numFrames, 120, 160, 1);
    lb_clip = max(lb, lb_min);
    ub_clip = min(ub, ub_max);

    % Create the volume star
    VS = VolumeStar(lb_clip, ub_clip);
end

function Xadv = PGD(net, X, T, varargin)
    p = inputParser;
    addParameter(p, 'Epsilon', 1/255);
    addParameter(p, 'Alpha', 2/255);
    addParameter(p, 'Steps', 10);
    addParameter(p, 'DataMin', 0);
    addParameter(p, 'DataMax', 1);

    parse(p, varargin{:});
    eps = p.Results.Epsilon;
    alpha = p.Results.Alpha;
    K = p.Results.Steps;
    xmin = p.Results.DataMin;
    xmax = p.Results.DataMax;

    X0 = X;

    % random start
    R = (2*rand(size(X0), 'like', extractdata(X0)) - 1) * eps;
    Xadv = X0 + dlarray(R, dimsLike(X0));
    Xadv = clipToRange(Xadv, xmin, xmax);
    Xadv = projectToLinfBall(Xadv, X0, eps);

    for k=1:K
        % compute gradient of loss wrt inputs
        [loss, gX] = dlfeval(@gradWrtInput, net, Xadv, T);
        
        step = alpha * sign(gX);
        Xadv = Xadv + step;

        Xadv = projectToLinfBall(Xadv, X0, eps);
        Xadv = clipToRange(Xadv, xmin, xmax);
    end
end

function [loss, gX] = gradWrtInput(dlnet, X, T)
    Y = forward(dlnet, X);
    loss = crossentropy(Y, T, 'TargetCategories', 'independent');
    gX = dlgradient(loss, X);
end

function Xc = clipToRange(X, a, b)
    Xc = min(max(X, a), b);
end

function Xp = projectToLinfBall(X, X0, eps)
    Xp = min(max(X, X0 - eps), X0 + eps);
end

function d = dimsLike(X)
    if hasdims(X)
        d = dims(X);
    else
        d = 'SSSCB';
    end
end