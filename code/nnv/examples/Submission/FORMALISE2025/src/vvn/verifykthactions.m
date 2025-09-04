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
    data = readNPY(sprintf("data/KTHActions/kthactions_%df.npy", smpLen));
    labels = readNPY(sprintf("data/KTHActions/kthactions_labels.npy"));

    % Preprocessing
    datacopy = data(:,:,:,:,:);

    % Experimental variables
    numClasses = 6;

    % Size of attack
    epsilon = [1/255; 2/255; 3/255];
    nE = length(epsilon);

    % Load the model
    modelName = sprintf("kthactions_%df.onnx", smpLen);
    netonnx = importONNXNetwork("models/" + modelName, "InputDataFormats", "TBCSS", "OutputDataFormats", "BC");
    net = matlab2nnv(netonnx);
    net.OutputSize = numClasses;
    disp("Finished loading model: " + modelName);

    % Verification settings
    reachOptions = struct;
    if verAlg == "relax"
        reachOptions.reachMethod = "relax-star-area";
        reachOptions.relaxFactor = 0.5;
    elseif verAlg == "approx"
        reachOptions.reachMethod = "approx-star";
    end

    % Make predictions on test set to check that we are verifying correct
    % samples
    outputLabels = zeros(length(datacopy));

    s = datacopy(index,:,:,:,:);
    s = squeeze(s);

    output = net.evaluate(s);
    [~, P] = max(output);

    %%%%%%%%%%%%%%%%
    % VERIFICATION %
    %%%%%%%%%%%%%%%%

    eps = epsilon(epsIndex);
    fprintf('Starting verification with epsilon %d \n', eps);
    
    % Get the sample
    sample = squeeze(datacopy(index,:,:,:,:));
    
    % Perform L_inf attack
    VS = L_inf_attack(sample, eps, smpLen);
    t = tic;

    % NEED THIS HERE SO MET EXISTS
    met = verAlg;
 
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
    lb_min = zeros(numFrames, 120, 160, 3);
    ub_max = ones(numFrames, 120, 160, 3);
    lb_clip = max(lb, lb_min);
    ub_clip = min(ub, ub_max);

    % Create the volume star
    VS = VolumeStar(lb_clip, ub_clip);
end

