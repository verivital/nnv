function [result,rT] = reach_acasxu(net, propertyFile)
%% Reachability analysis of acasxu benchmarks

    % Load specification to verify
    [lb_x, ub_x, property] = load_acasxu_vnnlib(propertyFile);
    % Create reachability parameters and options
    X = ImageStar(lb_x',ub_x');
    reachOptions = struct;
%     reachOptions.reachMethod = 'approx-star';
    reachOptions.reachMethod = 'approx-star';
    % Compute reachability
    rT = tic;
    Y = net.reach(X, reachOptions); % Seems to be working
    rT = toc(rT);
    % 
    result = verify_specification(Y, property); 
    % Evaluate property
    disp(' ');
    disp('===============================')
    disp('RESULTS')
    disp(' ')
    
    if result == 2
        if contains(net.reachMethod, "exact")
            result = 0;
            disp('Property is UNSAT');
        else
            disp('Property is UNKNOWN');
        end
    elseif result == 1
        disp('Property is SAT');
    else
        disp('Property is UNSAT')
    end
    
    disp("Reachability computation time = "+string(rT) + " seconds")

end