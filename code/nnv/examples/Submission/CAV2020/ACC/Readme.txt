
NOTE***: To run ACC case study, we need to have CORA installed. 

1) Install NNV: 

        a) git clone --recursive https://github.com/verivital/nnv.git

        b) go to ..code/nnv/ folder to run install.m

2) Run ACC case study: 

        a) Run ACC with discrete linear plant model (table 3 part 1): 

            run verify_linear_ACC.m 
            
            total runtime is about 90 seconds
            
        b) Run ACC with nonlinear plant model (table 3 part 2):

           run verify_nonlinear_ACC.m 

           total runtime is about 30 minutes

        c) To plot the reachable set figure (i.e. Figure 4)

           run plot_linear_ACC_reachSets.m 

           total runtime is about 60 seconds

NOTE***: 

We use 4 cores for computation. If your computer does not have >= 4 cores

Set numCores (or n_cores) parameter (default = 4) in the script to a suitable number 