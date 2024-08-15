function results = verify_instance_shape(net, vol, target, attack, reachOptions)
    % verify medmnist with inputs (input images), targets (labels) and attack
    % (struct with adversarial attack info)
    % results =  verify_medmnist(net, matlabNet, inputs, targets, attack, max_value*, min_value*)
    
    % Check what type of attack to consider
    if strcmp(attack.Name, 'add') || strcmp(attack.Name, 'remove')
        max_pixels = attack.max_pixels;
        noise_disturbance = attack.noise_de;
    else
        error("Adversarial attack not supported.");
    end

    % Choose attack
    if strcmp(attack.Name, 'add')
        I = add_voxels(vol, max_pixels, noise_disturbance);
    elseif strcmp(attack.Name, 'remove')
        I = remove_voxels(vol, max_pixels, noise_disturbance);
    end

    % Begin analysis

    t = tic; % start timer

    results = zeros(1,2);

    % Check for missclassification
    vol = single(vol);
    y = net.evaluate(vol);
    [~, y] = max(y);
    if y ~= target
        results(1) = -1; % missclassified
        results(2) = toc(t);
        return
    end

    % Check for falsification 
    n_samples = 100; % number of random samples to try for falsification
    xRand = I.sample(n_samples);
    for k = 1:n_samples
        x = xRand{k};
        y = net.evaluate(x);
        [~,idx] = max(y);
        if idx ~= target
            results(1) = 0;
            results(2) = toc(t);
            return
        end
    end

    % Compute reachability for verification
    try
        results(1) = net.verify_robustness(I, reachOptions, target);
    catch ME
        results(1) = -2;
        warning(ME.message);
    end

    % Save results
    results(2) = toc(t);

end


