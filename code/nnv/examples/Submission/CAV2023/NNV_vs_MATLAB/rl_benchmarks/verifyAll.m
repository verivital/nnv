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
    N_props = 50;
    res = zeros(N_props,6); % res, time, res, time, res, time
    % Specify NNV reachbility options
    reachOpt1 = struct; 
    reachOpt1.reachMethod = 'approx-star';
%     reachOpt2 = struct;
%     reachOpt2.reachMethod = 'exact-star';
    rng(1); % define random seed to reproduce results
    inst_idxs = randperm(200,N_props); % verify a subset of 296 properties

    % Begin reachability
    k = 1;
    for i=inst_idxs
        % Get network
        name = split(NNs_props_timeout.Var1{i},'/');
        name = name{2};
        net = networks{name2idx(name)};
        % Get vnnlib property
        propertyFile = string(NNs_props_timeout.Var2{i});
        % Run verification for each method
        [res(k,1), res(k,2)] = verify_rl_matlab(net.matlab, propertyFile);
        [res(k,3), res(k,4)] = verify_rl_nnv(net.nnv, propertyFile, reachOpt1);
%         [res(k,5), res(k,6)] = verify_rl_nnv(net.nnv, propertyFile, reachOpt2);
        k = k+1;
    end
    % Save results (including idxs of properties verified
    save("results_rl.mat","res", "inst_idxs");
end
