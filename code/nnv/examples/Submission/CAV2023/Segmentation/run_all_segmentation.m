function run_all_segmentation()
    
    %% run robustness analysis for MNIST networks

    cd MNIST;
    
    compare_mnist_nets_vs_num_attackedpixels; 
    compare_mnist_nets_vs_inputsize; % produce Figure 4a, Figure 5a, Figure 5b, Figure 8a
    
    compare_mnist_nets_vs_num_attackedpixels; % produce Figure 4b, Figure 8b
    
    compare_mnist_ReLU_reachTime_vs_others; % produce Figure 8c
    
    mnist_net1_verifyTime_vs_relaxFactor; % produce Table 2 - N1 results
    mnist_net2_verifyTime_vs_relaxFactor; % produce Table 2 - N2 results
    mnist_net3_verifyTime_vs_relaxFactor; % produce Table 2 - N3 results
        
    %% run robustness analysis for M2NIST networks
    
    cd M2NIST;

    compare_m2nist_nets_vs_num_attackedpixels; % produce Figure 6a, Figure 6b, Figure 7a, Figure 7b
    

end