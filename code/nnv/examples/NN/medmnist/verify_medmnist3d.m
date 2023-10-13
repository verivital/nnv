function results = verify_medmnist3d(net, matlabNet, inputs, targets, attack, max_value, min_value)
    % verify medmnist with inputs (input images), targets (labels) and attack
    % (struct with adversarial attack info)
    % results =  verify_medmnist(net, matlabNet, inputs, targets, attack, max_value*, min_value*)
    
    % Check what type of attack to consider
    if strcmp(attack.Name, 'linf')
        epsilon = attack.epsilon;
        max_pixels = attack.max_pixels;
    else
        error("Adversarial attack not supported.");
    end

    % check if max_value and min_value are provided
    if ~exist("max_value", 'var')
        max_value = inf;
    end
    if ~exist("min_value", 'var')
        min_value = -inf;
    end

    % check the inputs size is in the correct order (num of volumes in dim 4)
    N = length(targets); % number of volumes to verify
    if size(inputs,5) ~= N
        error("Number of volumes and targets must be the same. Ensure number of volumes are in dimension 5 of input array.")
    end

    % Initialize variables
    results = zeros(2,N); % 1-> verify result, 2-> verify time
    
    % Define reachability parameters
    reachOptions = struct;
    reachOptions.reachMethod = 'approx-star';
    reachOptions.dis_opt = 'display';
    
    % Evaluate all images
    for i = 1:N

        % print progress
        % if ~mod(i, 3)
        disp("Verifying input volume " + string(i) + " of " + string(N) + "...");
        % end

        % Create set of volume images
        vol = inputs(:,:,:,:,i);
        I = L_inf_attack(vol, epsilon, max_pixels, max_value, min_value);

        t = tic; % start timer

        % Check for missclassification
        vol = single(vol);
        y = net.evaluate(vol);
        [~, y] = max(y);
        if y ~= targets(i)
            results(1, i) = -1; % missclassified
            results(2,i) = toc(t);
            continue;
        end

        % Check for falsification with upper and lower bounds
        yUpper = net.evaluate(I.vol_ub);
        [~, yUpper] = max(yUpper);
        yLower = net.evaluate(I.vol_lb);
        [~, yLower] = max(yLower);
        if yUpper ~= targets(i) || yLower ~= targets(i)
            results(1,i) = 0; % not robust
            results(2,i) = toc(t);
            continue;
        end

        % Compute reachability for verification
        results(1,i) = net.verify_robustness(I, reachOptions, targets(i));
        results(2,i) = toc(t);

    end

end

% Return a VolumeStar of an Linf attack
function I = L_inf_attack(vol, epsilon, max_pixels, max_value, min_value)
    volSize = size(vol);
    n = numel(vol);
    idxs = randperm(n, max_pixels); % select some random pixels from the input volume
    disturbance = zeros(volSize, "like", vol); % initialize array (mostly zeros)
    disturbance(idxs) = epsilon; % change to {epsilon} the randomly selected pixels
    lb = max(vol - disturbance, min_value);
    ub = min(vol + disturbance, max_value);
    I = VolumeStar(single(lb), single(ub)); % default: single (assume onnx input models)
end

