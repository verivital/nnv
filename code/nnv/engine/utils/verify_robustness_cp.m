function result = verify_robustness_cp(net,inputSet, reachOptions, target)
    % net: dlnetwork
    % Compute reachable set
    R = Prob_reach(net, inputSet, reachOptions);
    % Check robustness
    result = obj.checkRobust(R, target);
end

