function results = verify_medmnist2d(net, inputs, targets, attack, max_value, min_value)
    % verify medmnist with inputs (input images), targets (labels) and attack
    % (struct with adversarial attack info)
    % results =  verify_medmnist(inputs, targets, attack, max_value*, min_value*)
    
    % Check what type of attack to consider
    if strcmp(attack.Name, 'linf')
        epsilon = attack.epsilon;
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

    % check the inputs size is in the correct order (num of images in dim 4)
    N = length(targets); % number of images to verify
    if size(inputs,4) ~= N
        error("Number of images and targets must be the same. Ensure number of images are in dimension 4 of inputs.")
    end

    % Initialize variables
    results = zeros(2,N); % 1-> verify result, 2-> verify time
    
    % Define reachability parameters
    reachOptions = struct;
    reachOptions.reachMethod = 'approx-star';
    
    % Evaluate all images
    for i = 1:N

        % print progress
        % if ~mod(i, 3)
        disp("Verifying image " + string(i) + " of " + string(N) + "...");
        % end

        % Create set of images
        img = inputs(:,:,:,i);
        I = l_inf_attack(img, epsilon, max_value, min_value);

        t = tic; % start timer

        % Check for missclassification
        img = single(img);
        y = net.evaluate(img);
        [~, y] = max(y);
        if y ~= targets(i)
            results(1, i) = -1; % missclassified
            results(2,i) = toc(t);
            continue;
        end

        % Check for falsification with upper and lower bounds
        yUpper = net.evaluate(I.im_ub);
        [~, yUpper] = max(yUpper);
        yLower = net.evaluate(I.im_lb);
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

% Return an ImageStar of an linf attack
function I = l_inf_attack(img, epsilon, max_value, min_value)
    imgSize = size(img);
    disturbance = epsilon * ones(imgSize, "like", img); % disturbance value
    lb = max(img - disturbance, min_value);
    ub = min(img + disturbance, max_value);
    I = ImageStar(single(lb), single(ub)); % default: single (assume onnx input models)
end

