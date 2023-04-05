function reach_random(neuralode,x0,fname,unc)

    % Create initial state
    lb = x0-unc;
    ub = x0+unc;
    R0 = Star(lb,ub);
    
    % y = neuralode.evaluate(x0); % Simulation
    t = tic;
    R = neuralode.reach(R0); % Reachability
    rT = toc(t);
    
    % save results
    rname = string(fname)+"_"+string(unc)+".mat";
    save(rname,'R','rT');

end

