function verifyAll()
    
    % Load network/property combinations
    csvFile = "instances.csv";
    opts = detectImportOptions(csvFile);
    opts.Delimiter = ',';
    NNs_props_timeout = readtable(csvFile, opts);
    N = height(NNs_props_timeout);
    % Load networks
    [networks, name2idx] = load_rl_NNs();
    % Init result variable
    res = zeros(N,6); % res, time, res, time, res, time
    % Specify NNV reachbility options
    reachOpt1 = struct; 
    reachOpt1.reachMethod = 'approx-star';
    reachOpt2 = struct;
    reachOpt2.reachMethod = 'exact-star';
    max_cores = getenv('NUMBER_OF_PROCESSORS');
    reachOpt2.numCores = min(8,max_cores); % try 8, but if local system has less, compute with max number of cores

    % Begin reachability
    for i=1:N
        % Get network
        name = split(NNs_props_timeout.Var1{i},'/');
        name = name{2};
        net = networks{name2idx(name)};
        % Get vnnlib property
        propertyFile = string(NNs_props_timeout.Var2{i});
        % Run verification for each method
        [res(i,1), res(i,2)] = verify_rl_matlab(net.matlab, propertyFile);
        [res(i,3), res(i,4)] = verify_rl_nnv(net.nnv, propertyFile, reachOpt1);
        [res(i,5), res(i,6)] = verify_rl_nnv(net.nnv, propertyFile, reachOpt2);
    end
    % Save results
    save("results_rl.mat","res");
end
