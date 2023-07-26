
1. install the NNV tool and Add path to this directory

2. Download datasets: 
    i. Battery State-of-Charge Dataset(BSOC) : follow 'BatterySOCUsingDeepLearningExample.mlx' file for downloading
                                        datasets and necessary pre-processing
    ii. Turbine Engine Degradation Simulation (TEDS) : follow 'RULEstimationUsingCNNExample' file for downloading
                                        datasets and necessary pre-processing

3. File/Folder details:
     a. Appendix : Contains all the images shown in the Appendix Section
            i. BatterySOC dataset   : Fig.7 in paper, sample dataset plot for BSOC dataset
           ii. TEDS_RULdataset      : Fig.8 in paper, sample dataset plot for TEDS dataset
          iii. Local_Robustness11   : Fig.6 in paper, example for calculating different robustness measure 
           iv. Network_Architecture : Fig.9 in paper, the networks used for the two datasets
            v. noise_MFAI           : Fig.10 in paper, spread of noise in case of MFAI noise
           vi. noise_MFSI           : Fig.11 in paper, spread of noise in case of MFSI noise
          vii. noise_SFAI           : Fig.12 in paper, spread of noise in case of SFAI noise
         viii. noise_SFSI           : Fig.13 in paper, spread of noise in case of SFSI noise
           ix. RULPlot_SFAI         : Fig.14 in paper, reachability plot for TEDS dataset for a particular time-segment with increasing noise

     b. PredictBatterySOCUsingDeepLearningExample : Main folder for running the tests for BSOC dataset
            i. BatterySOCUsingDeepLearningExample : for data download, pre-process and network training
           ii. createSOCDataset                   : function for creating time series dataset for reachability analysis
          iii. reachabilityNoiseForSOC            : function for adding pertubations in the input time series data
           iv. generate_paper_results             : contains commands for generating results shown in the paper.
            v. SOC_results                        : all the result variables stores in a .mat file
           vi. generate_fig3_robustness_value_SOC : generates plot for robustness measures with increasing noise (Fig.3)
          vii. generate_fig2_BSOC_SFAI            : reachability plot for BSOC dataset for a particular time-segment with increasing noise (Fig.2)

     c. RULEstimationUsingCNNExample              : Main folder for running the tests for TEDS dataset
            i. RULEstimationUsingCNNExample : for data download, pre-process and network training
           ii. createTEDSDataset            : function for creating time series dataset for reachability analysis
          iii. reachabilityNoiseForTEDS     : function for adding pertubations in the input time series data
           iv. generate_paper_results       : contains commands for generating results shown in the paper.
            v. TEDS_results                 : all the result variables stores in a .mat file
           vi. plot_localmonotonicity       :
          vii. localmonotonicity            : plot for local monotonicity shown in paper (Fig.4)
         viii. generate_fig5_robustness_value_RUL: generates plot for robustness measures with increasing noise (Fig.5)
    
    d. localrobustnessSingle            : function for calculating robustness value at a time-instance
    e. localrobustnessSingle_overlap    : function for calculating the overlap between the estimated value and the permissible value at a time-instance
    f. globalrobustnessSingle           : function for calculating the percentage sample robustness
    g. globalrobustnessSingle_overlap   : function for calculating percentage overlap robustness
    h. plot_robustnessvalues            :

    i. Star         : the star image shown in paper (Fig.1)
    