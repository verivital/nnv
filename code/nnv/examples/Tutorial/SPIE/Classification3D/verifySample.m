function [res] = verifySample(net, I, img, target, reachOptions)
    
    % Check for missclassification
    img = single(img);
    y = net.evaluate(img);
    [~, y] = max(y);
    if y ~= target
        res = -1; % missclassified
        return;
    end

    % Check for falsification with upper and lower bounds
    yUpper = net.evaluate(I.vol_ub);
    [~, yUpper] = max(yUpper);
    yLower = net.evaluate(I.vol_lb);
    [~, yLower] = max(yLower);
    if yUpper ~= target || yLower ~= target
        res = 0; % not robust
        return;
    end

    % Compute reachability for verification
    try
        res = net.verify_robustness(I, reachOptions, target);
    catch ME
        warning(ME.message);
        res = -2;
    end

end

