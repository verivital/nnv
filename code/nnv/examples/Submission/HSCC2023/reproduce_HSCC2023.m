% produce all results of HSCC2023

pwd;
cd small_RNN;

%% run robustness verification of the small RNN network
verify_small_RNN; % produce Figure 3, Table 1, Figure 4.


cd ../N_2_0;
%% run robustness verification of the N_2_0 RNN network
clear;
verify_N2_0; % produce N_2_0 part of Table 2 
clear;
verify_N2_0_full; % produce Table 4


cd ../N_2_2;
%% run robustness verification of the N_2_2 RNN network
clear;
verify_N_2_2; % produce N_2_2 part of Table 2 
clear;
verify_N_2_2_full; % produce Table 5


cd ../N_4_0;
%% run robustness verification of the N_4_0 RNN network
clear;
verify_N4_0; % produce N_4_0 part of Table 2 
clear;
verify_N4_0_full; % produce Table 6


cd ../N_4_2;
%% run robustness verification of the N_4_2 RNN network
clear;
verify_N4_2; % produce N_4_2 part of Table 2 
clear;
verify_N4_2_full; % produce Table 7


cd ../N_4_4;
%% run robustness verification of the N_4_4 RNN network
clear;
verify_N4_4; % produce N_4_4 part of Table 2 
clear;
verify_N4_4_full; % produce Table 8



cd ../N_8_0;
%% run robustness verification of the N_8_0 RNN network
clear;
verify_N8_0; % produce N_8_0 part of Table 2 
clear;
verify_N8_0_full; % produce Table 9
clear;
epsilon_vs_robustness_and_time; % produce Figure 5
