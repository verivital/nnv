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
    res = zeros(N,2); % res, time
    % Specify NNV reachbility options
    reachOpt1 = struct; 
    reachOpt1.reachMethod = 'approx-star';

    % Begin reachability
    for i=1:N
        % Get network
        name = split(NNs_props_timeout.Var1{i},'/');
        name = name{2};
        net = networks{name2idx(name)};
        % Get vnnlib property
        propertyFile = string(NNs_props_timeout.Var2{i});
        % Run verification for each method
        [res(i,1), res(i,2)] = verify_oval21_nnv(net, propertyFile, reachOpt1);
    end
    
    % Save results
    save("results_oval21.mat","res");

end
