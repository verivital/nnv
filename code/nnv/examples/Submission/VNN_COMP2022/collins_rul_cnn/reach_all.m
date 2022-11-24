%% Run all properties of collins benchmarks
% benchmarkFolder = "/home/manzand/Documents/MATLAB/vnncomp2022_benchmarks/benchmarks/collins_rul_cnn/";
benchmarkFolder = "/home/dieman95/Documents/MATLAB/vnncomp2022_benchmarks/benchmarks/collins_rul_cnn/";

%% 1) Load the networks
[networks, name2idx] = load_collins_NNs();

%% 2) Verify all properties
csvFile = "instances.csv";
NNs_props_timeout = readtable(benchmarkFolder+csvFile);
% run each instance in MATLAB
% ideally, for the competition, we'll create a matlab function that
% iterates through each csv and sets a timeout for each instance (column 3)
verified = 0;
verfied = zeros(height(NNs_props_timeout),1);
errorMsgs = [];
reachTimes = zeros(height(NNs_props_timeout),1);
for i=1:height(NNs_props_timeout)
    name = split(NNs_props_timeout.Var1{i},'/');
    name = name{2};
    net = networks{name2idx(name)};
    propertyFile = benchmarkFolder + string(NNs_props_timeout.Var2{i});
    try
        [result,rT] = reach_collins(net,propertyFile);
        verified(i) = result;
        reachTimes(i) = rT;
    catch ME
        verified(i) = -1; % verification failed
        errorMsgs = [errorMsgs; ME];
        warning("Reachability of network " + string(name) + " WITH specification " + string(NNs_props_timeout.Var2{i}) + " FAILED");
    end
end

%% Show results
disp(' ')
disp('=============================')
disp('collins_rul_cnn Benchmarks');
disp('=============================')
disp('Total number of properties verified = '+string(sum(verified==1)) + " out of "+height(NNs_props_timeout));
disp("UNSAT = "+string(sum(verified==0)))
disp("Error/timeout = "+string(sum(verified==-1)))
disp("Total computation time = " + string(sum(reachTimes)) + " seconds")
% Could take a look into the unsat ones and redo the reachability analysis
% with exact star