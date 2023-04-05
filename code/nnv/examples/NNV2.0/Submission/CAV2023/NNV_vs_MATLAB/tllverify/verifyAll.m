function verifyAll()
    
    % Load network/property combinations
    csvFile = "instances.csv";
    opts = detectImportOptions(csvFile);
    opts.Delimiter = ',';
    NNs_props_timeout = readtable(csvFile, opts);
    N = 10;
    % Init result variable
    res = zeros(N,6); % res, time, res, time, res, time
    
    % Specify NNV reachbility options
    reachOpt1 = struct; 
    reachOpt1.reachMethod = 'approx-star';

    % Reachability computation
    for i=1:N
        [res(i,1), res(i,2)] = verify_tllverify_matlab(NNs_props_timeout.Var1{i}, NNs_props_timeout.Var2{i});         % MATLAB
        [res(i,3), res(i,4)] = verify_tllverify_nnv(NNs_props_timeout.Var1{i}, NNs_props_timeout.Var2{i}, reachOpt1); % approx-star
    end

    % Save results
    save('results_tllverify.mat', 'res');

end