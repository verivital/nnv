function verifyAll()
    
    % Load network/property combinations
    csvFile = "instances.csv";
    opts = detectImportOptions(csvFile);
    opts.Delimiter = ',';
    NNs_props_timeout = readtable(csvFile, opts);
    N = height(NNs_props_timeout);
    % Load networks
    [networks, name2idx] = load_oval21_NNs();
    % Init result variable
    res = zeros(N,4); % res, time, res, time
    % Specify NNV reachbility options
    reachOpt1 = struct; 
    reachOpt1.reachMethod = 'approx-star';
%     reachOpt2 = struct;
%     reachOpt2.reachMethod = 'exact-star';
%     reachOpt2.lp_solver = 'glpk';
%     max_cores = getenv('NUMBER_OF_PROCESSORS');
%     reachOpt2.numCores = min(8,max_cores); % try 8, but if local system has less, compute with max number of cores
    rng(0); % define random seed to reproduce results
%     inst_idxs = randperm(N,2); % verify just 2 properties out of 30
    inst_idxs = [15,19];

    % Begin reachability
    for i=inst_idxs
        % Get network
        name = split(NNs_props_timeout.Var1{i},'/');
        name = name{2};
        net = networks{name2idx(name)};
        % Get vnnlib property
        propertyFile = string(NNs_props_timeout.Var2{i});
        % Run verification for each method
        [res(i,1), res(i,2)] = verify_oval21_nnv(net, propertyFile, reachOpt1);
%         [res(i,3), res(i,4)] = verify_oval21_nnv(net, propertyFile, reachOpt2);
    end
    % Save results
    save("results_oval21.mat","res");
end
