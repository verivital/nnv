%% Process acas xu results
% Load results
p1 = load("results_p1.mat");
p2 = load("results_p2.mat");
p3 = load("results_p3.mat");
p4 = load("results_p4.mat");

% Process them
% Property 1
p1_ver = sum(p1.resNNV == categorical("verified"));
p1_vio = sum(p1.resNNV == categorical("violated"));
p1_unk = sum(p1.resNNV == categorical("unproven"));
% Property 2
p2_ver = sum(p2.resNNV == categorical("verified"));
p2_vio = sum(p2.resNNV == categorical("violated"));
p2_unk = sum(p2.resNNV == categorical("unproven"));
% Property 3
p3_ver = sum(p3.resNNV == categorical("verified"));
p3_vio = sum(p3.resNNV == categorical("violated"));
p3_unk = sum(p3.resNNV == categorical("unproven"));
% Property 4
p4_ver = sum(p4.resNNV == categorical("verified"));
p4_vio = sum(p4.resNNV == categorical("violated"));
p4_unk = sum(p4.resNNV == categorical("unproven"));
% Results do not match those of other tools, take a look at how we are
% verifying the specifications