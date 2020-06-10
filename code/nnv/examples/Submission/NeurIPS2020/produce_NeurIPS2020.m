% produce all results of NeurIPS2020 

pwd; 
cd MNIST;

% run robustness analysis for MNIST networks
compare_mnist_nets_vs_num_attackedpixels; 
clear;
compare_mnist_nets_vs_inputsize;


mydir  = pwd;
idcs   = strfind(mydir,filesep);
newdir = mydir(1:idcs(end)-1);
cd M2NIST;

% run robustness analysis for M2NIST networks
clear;
compare_m2nist_nets_vs_num_attackedpixels;