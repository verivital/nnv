%% Run all properties of acas xu benchmarks
% benchmarkFolder = "/home/manzand/Documents/MATLAB/vnncomp2022_benchmarks/benchmarks/acasxu/";
benchmarkFolder = "/home/dieman95/Documents/MATLAB/vnncomp2022_benchmarks/benchmarks/acasxu/";

%% 1) Load the networks
% [networks, name2idx] = load_acasxu_NNs();

%% 2) Verify all properties
csvFile = "instances.csv";
opts = detectImportOptions(benchmarkFolder+csvFile);
opts.Delimiter = ',';
NNs_props_timeout = readtable(benchmarkFolder+csvFile, opts);
% run each instance in MATLAB
% ideally, for the competition, we'll create a matlab function that
% iterates through each csv and sets a timeout for each instance (column 3)
verified = 0;
verfied = zeros(height(NNs_props_timeout),1);
errorMsgs = [];
reachTimes = zeros(height(NNs_props_timeout),1);
for i=1:height(NNs_props_timeout)
    name = split(NNs_props_timeout.Var1{i}, filesep);
    name = name{2};
    net = networks{name2idx(name)};
    propertyFile = benchmarkFolder + string(NNs_props_timeout.Var2{i});
    try
        [result,rT] = reach_acasxu(net,propertyFile);
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
disp('ACAS Xu Benchmarks');
disp('=============================')
disp('Total number of properties = ' + string(height(NNs_props_timeout)));
disp("SAT = " +string(sum(verified==1)));
disp("UNSAT = "+string(sum(verified==0)));
disp("UNKNOWN = " + string(sum(verified==2)));
disp("Error/timeout = "+string(sum(verified==-1)));
disp("Total computation time = " + string(sum(reachTimes)) + " seconds")
% Could take a look into the unsat ones and redo the reachability analysis
% with exact star