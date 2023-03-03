% produce all results of NeurIPS2021 

pwd; 
cd MNIST;

%% run robustness analysis for MNIST networks
compare_mnist_nets_vs_num_attackedpixels; 
clear;
compare_mnist_nets_vs_inputsize; % produce Figure 4a, Figure 5a, Figure 5b, Figure 8a

clear; 
compare_mnist_nets_vs_num_attackedpixels; % produce Figure 4b, Figure 8b

clear;
compare_mnist_ReLU_reachTime_vs_others; % produce Figure 8c

clear;
mnist_net1_verifyTime_vs_relaxFactor; % produce Table 2 - N1 results
clear;
mnist_net2_verifyTime_vs_relaxFactor; % produce Table 2 - N2 results
clear;
mnist_net3_verifyTime_vs_relaxFactor; % produce Table 2 - N3 results


mydir  = pwd;
idcs   = strfind(mydir,filesep);
newdir = mydir(1:idcs(end)-1);
cd M2NIST;

%% run robustness analysis for M2NIST networks
clear;
compare_m2nist_nets_vs_num_attackedpixels; % produce Figure 6a, Figure 6b, Figure 7a, Figure 7b


mydir  = pwd;
idcs   = strfind(mydir,filesep);
newdir = mydir(1:idcs(end)-1);
cd RelaxStarCompare;

%% run comparison of different relaxation strategies 
clear;
mnist01_conservativeness_vs_relaxFactor; % produce Figure 9 (combined from 2 figures)

clear;
mnist01_relax_star_area_vs_eran; % produce Figure 10 (combined from 2 figures)
