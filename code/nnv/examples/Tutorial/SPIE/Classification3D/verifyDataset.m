function results = verifyDataset(net, inputs, targets, attack, max_value, min_value)
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
    % reachOptions.reachMethod = 'approx-star';
    reachOptions.reachMethod = 'relax-star-area';
    reachOptions.relaxFactor = 1;
    
    % Analyze all images
    for i = 1:N

        % print progress
        % if ~mod(i, 3)
        disp("Verifying image " + string(i) + " of " + string(N) + "...");
        % end

        % Create set of images
        img = inputs(:,:,:,i);
        I = l_inf_set(img, epsilon, max_value, min_value);
        target = targets(i);

        t = tic; % start timer
        results(1,i) = verifySample(net,I,img,target,reachOptions);
        results(2,i) = toc(t);

    end

end



