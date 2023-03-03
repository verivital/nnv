if not(isfolder('results'))
    mkdir results; % Create directory to store results
end
t_total = tic;

% Run ACC 
cd benchmarks/ACC;
if not(isfolder('../../results/ACC'))
    mkdir ../../results/ACC; % Create directory to store results
end
try
    run reach.m
catch e %e is an MException struct
    fprintf(1,'The identifier was:\n%s',e.identifier);
    fprintf(1,'There was an error! The message was:\n%s',e.message);
end

% Run Airplane
cd ../Airplane;
if not(isfolder('../../results/Airplane'))
    mkdir ../../results/Airplane; % Create directory to store results
end
try
    run reach.m
catch e %e is an MException struct
    fprintf(1,'The identifier was:\n%s',e.identifier);
    fprintf(1,'There was an error! The message was:\n%s',e.message);
end

% Run Benchmark 9 - Tora
cd ../Benchmark9-Tora;
if not(isfolder('../../results/benchmark9-tora'))
    mkdir ../../results/benchmark9-tora; % Create directory to store results
end
try
    run reach.m
catch e %e is an MException struct
    fprintf(1,'The identifier was:\n%s',e.identifier);
    fprintf(1,'There was an error! The message was:\n%s',e.message);
end

% Run Benchmark10 = Unicycle
cd ../Benchmark10-Unicycle;
if not(isfolder('../../results/Unicycle'))
    mkdir ../../results/Unicycle; % Create directory to store results
end
try
    run reach.m
catch e %e is an MException struct
    fprintf(1,'The identifier was:\n%s',e.identifier);
    fprintf(1,'There was an error! The message was:\n%s',e.message);
end

% Run Double Pendulum
cd ../Double_Pendulum;
if not(isfolder('../../results/DoublePendulum'))
    mkdir ../../results/DoublePendulum; % Create directory to store results
end
try
    run reach_dp_less.m;
catch e %e is an MException struct
    fprintf(1,'The identifier was:\n%s',e.identifier);
    fprintf(1,'There was an error! The message was:\n%s',e.message);
end
try
    run reach_dp_more.m;
catch e %e is an MException struct
    fprintf(1,'The identifier was:\n%s',e.identifier);
    fprintf(1,'There was an error! The message was:\n%s',e.message);
end

% Run Single Pendulum
cd ../Single_Pendulum;
if not(isfolder('../../results/SinglePendulum'))
    mkdir ../../results/SinglePendulum; % Create directory to store results
end
try
    run reach.m
catch e %e is an MException struct
    fprintf(1,'The identifier was:\n%s',e.identifier);
    fprintf(1,'There was an error! The message was:\n%s',e.message);
end

% Run Tora Heterogeneous
cd ../Tora_Heterogeneous;
if not(isfolder('../../results/Tora_Heterogeneous'))
    mkdir ../../results/Tora_Heterogeneous; % Create directory to store results
end
try
    run reachTora_reluTanh.m;
catch e %e is an MException struct
    fprintf(1,'The identifier was:\n%s',e.identifier);
    fprintf(1,'There was an error! The message was:\n%s',e.message);
end
try
    run reachTora_sigmoid.m;
catch e %e is an MException struct
    fprintf(1,'The identifier was:\n%s',e.identifier);
    fprintf(1,'There was an error! The message was:\n%s',e.message);
end

%% Run VCAS
cd ../VCAS;
if not(isfolder('../../results/VCAS'))
    mkdir ../../results/VCAS; % Create directory to store results
end

path_out = [path_results(), filesep, 'VCAS', filesep];
mkdir(path_out)

run reachVCAS_middle19.m
run reachVCAS_middle22.m
run reachVCAS_middle25.m
run reachVCAS_middle28.m
run reachVCAS_worst19.m
run reachVCAS_worst22.m
run reachVCAS_worst25.m
run reachVCAS_worst28.m

toc(t_total); % Time to run all experiments
cd ../..;