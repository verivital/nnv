function processResults()

%% Process acas xu results

    % Load results
    p3 = load("results_p3.mat");
    p4 = load("results_p4.mat");
    p3mat = load("results_p3_mat.mat");
    p4mat = load("results_p4_mat.mat");
    
    % Property 3
    % NNV
    p3_unsat_nnv = sum(p3.resNNV == categorical("verified"));
    p3_sat_nnv   = sum(p3.resNNV == categorical("violated"));
    p3_unk_nnv   = sum(p3.resNNV == categorical("nope"));
    p3_time_nnv  = sum(p3.timeNNV)/45;
    % MATLAB
    p3_sat_mat   = sum(p3mat.resMAT == categorical("verified"));
    p3_unsat_mat = sum(p3mat.resMAT == categorical("violated"));
    p3_unk_mat   = sum(p3mat.resMAT == categorical("unproven"));
    p3_time_mat  = sum(p3mat.timeMAT)/45;
    
    % Property 4
    p4_unsat_nnv = sum(p4.resNNV == categorical("verified"));
    p4_sat_nnv   = sum(p4.resNNV == categorical("violated"));
    p4_unk_nnv   = sum(p4.resNNV == categorical("nope"));
    p4_time_nnv = sum(p4.timeNNV)/45;
    % MATLAB
    p4_sat_mat   = sum(p4mat.resMAT == categorical("verified"));
    p4_unsat_mat = sum(p4mat.resMAT == categorical("violated"));
    p4_unk_mat   = sum(p4mat.resMAT == categorical("unproven"));
    p4_time_mat  = sum(p4mat.timeMAT)/45;

    % Create table from here

end