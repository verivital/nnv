function verifyAll()
    % Load network/property combinations
    csvFile = "instances.csv";
    opts = detectImportOptions(csvFile);
    opts.Delimiter = ',';
    NNs_props_timeout = readtable(csvFile, opts);
%     N = height(NNs_props_timeout);
    N = 4;
    % Init result variable
    res = zeros(N,6); % res, time, res, time, res, time
    % Specify NNV reachbility options
    reachOpt1 = struct; 
    reachOpt1.reachMethod = 'approx-star';
    reachOpt2 = struct;
    reachOpt2.reachMethod = 'exact-star';
    max_cores = getenv('NUMBER_OF_PROCESSORS');
    reachOpt2.numCores = min(8,max_cores); % try 8, but if local system has less, compute with max number of cores
    for i=1:N
        [res(i,1), res(i,2)] = verify_tllverify_matlab(NNs_props_timeout.Var1{i}, NNs_props_timeout.Var2{i});         % MATLAB
        [res(i,3), res(i,4)] = verify_tllverify_nnv(NNs_props_timeout.Var1{i}, NNs_props_timeout.Var2{i}, reachOpt1); % approx-star
        [res(i,5), res(i,6)] = verify_tllverify_nnv(NNs_props_timeout.Var1{i}, NNs_props_timeout.Var2{i}, reachOpt2); % exact-star
    end
    save('results_tllverify.mat', 'res');
end

% Example to load networks
%     net = importONNXNetwork("onnx/tllBench_n=2_N=M=8_m=1_instance_0_0.onnx", InputDataFormats="BC");
%     loadOpt.InputDataFormat = "BC";
%     nn = onnx2nnv("onnx/tllBench_n=2_N=M=8_m=1_instance_0_0.onnx", loadOpt);