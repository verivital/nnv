function paper_results(resVersion)

    if ~exist("resVersion", "var") % default is the short version
        resVersion = "short";
    end

    if strcmp(resVersion, "short")
        pShort;
    elseif strcmp(resVersion, "all")
        pLong;
    else % assume it's the long version
        error("Wrong resVersion value ({short, all}");
    end

end

%% Helper functions

% Create ACAS Xu comparison table
function acas_compare()
    % ACAS Xu comparison
    max_cores = maxNumCompThreads;
    used_cores = min([8,max_cores]);
    % Property 3

    % NNV
    p3 = load("NNV_vs_MATLAB/acas/results_p3.mat");
    p4 = load("NNV_vs_MATLAB/acas/results_p4.mat");
    p3_unsat_nnv = sum(p3.resNNV == categorical("verified"));
    p3_sat_nnv   = sum(p3.resNNV == categorical("violated"));
    p3_unk_nnv   = sum(p3.resNNV == categorical("nope"));
    p3_time_nnv  = sum(p3.timeNNV)/45;
    % MATLAB
    p3_sat_mat   = sum(p3.resMAT == categorical("verified"));
    p3_unsat_mat = sum(p3.resMAT == categorical("violated"));
    p3_unk_mat   = sum(p3.resMAT == categorical("unproven"));
    p3_time_mat  = sum(p3.timeMAT)/45;
    
    % Property 4
    p4_unsat_nnv = sum(p4.resNNV == categorical("verified"));
    p4_sat_nnv   = sum(p4.resNNV == categorical("violated"));
    p4_unk_nnv   = sum(p4.resNNV == categorical("nope"));
    p4_time_nnv = sum(p4.timeNNV)/45;
    % MATLAB
    p4_sat_mat   = sum(p4.resMAT == categorical("verified"));
    p4_unsat_mat = sum(p4.resMAT == categorical("violated"));
    p4_unk_mat   = sum(p4.resMAT == categorical("unproven"));
    p4_time_mat  = sum(p4.timeMAT)/45;

    % Now create the tables
    % TABLE 2 - ACAS Xu
    if is_codeocean()
        fid = fopen([path_results_codeocean, 'Table_2.txt'],'w'); 
    else
        fid = fopen('Table_2.txt','w'); 
    end
    fprintf(fid, '-------------------------------------------------------------------------------------------------\n');
    fprintf(fid, "Property 3 (45) \n");
    fprintf(fid, '-------------------------------------------------------------------------------------------------\n');
    fprintf(fid, "            matlab    approx    relax 25%%    relax 50%%    relax 75%%    relax 100%%    exact (%d)\n",used_cores);
    fprintf(fid, " SAT          %d         %d           %d            %d            %d             %d           %d \n", p3_sat_mat, p3_sat_nnv);
    fprintf(fid, " UNSAT        %d        %d          %d            %d            %d             %d           %d \n", p3_unsat_mat, p3_unsat_nnv);
    fprintf(fid, " time (s)   %.4g    %.4g     %.4g        %.4g       %.4g       %.4g       %.4g \n\n",p3_time_mat, p3_time_nnv );
    fprintf(fid, '-------------------------------------------------------------------------------------------------\n');
    fprintf(fid, "Property 4 (45) \n");
    fprintf(fid, '-------------------------------------------------------------------------------------------------\n');
    fprintf(fid, "            matlab    approx    relax 25%%    relax 50%%    relax 75%%    relax 100%%    exact (%d)\n",used_cores);
    fprintf(fid, " SAT          %d         %d           %d            %d            %d             %d           %d \n", p4_sat_mat, p4_sat_nnv);
    fprintf(fid, " UNSAT        %d         %d          %d            %d            %d             %d           %d \n", p4_unsat_mat, p4_unsat_nnv);
    fprintf(fid, " time (s)   %.4g    %.4g     %.4g        %.4g       %.4g       %.4g       %.4g \n", p4_time_mat, p4_time_nnv );
    fclose(fid);
end

% Compute total robustness of rnn
function rob = compute_robustness(rb)
    % Compute RNN robutness
    [n,m] = size(rb);
    rob = 0;
    for i=1:n
        for j=1:m
            if rb{i,j}(end) == 1
                rob = rob+1;
            end
        end
    end
end

% Short RE package
function pShort()

    % Results from section 4.2, 4.2, 4.4

    %FPA
    rACC = load("NeuralODEs/ACC/results_tanhplant.mat");
    % FPA
    rFPA = load("NeuralODEs/FPA/fpa_reach.mat");
    % MNIST
    rMid = load("NeuralODEs/MNIST/cnn_medium_nnv_inf_0.5.mat");
    rTiny = load("NeuralODEs/MNIST/cnn_tiny_nnv_inf_0.5.mat");
    
    % RNN
    rR2 = load("RNN/N_2_0/N2_0_results.mat");
    rR4 = load("RNN/N_4_4/N4_4_results.mat");
    rR8 = load("RNN/N_8_0/N8_0_results.mat");

    % Segmentation
    rDil = load("Segmentation/dilated_results_0.0001.mat");
    rTrans = load("Segmentation/transposed_results_0.0001.mat");
    
    % Write all results to a text file
    if is_codeocean()
        fid = fopen([path_results_codeocean, 'results_4.2-4.4.txt'],'w'); 
    else
        fid = fopen('results_4.2-4.4.txt','w'); 
    end

    % Section 4.2
    fprintf(fid, "-------------------------------------\n");
    fprintf(fid, "SECTION 4.2 - Neural ODEs \n");
    fprintf(fid, "-------------------------------------\n\n");
    fprintf(fid, "Dynamical system - FPA\n");
    fprintf(fid, "Computation time: %f seconds.\n\n", rFPA.time);
    fprintf(fid, "Classification - MNIST\n");
    fprintf(fid, "CNODE_{S} robustness: %f with average computation time: %f seconds.\n", sum(rTiny.rob_ode)/rTiny.numT, rTiny.timeT/rTiny.numT);
    fprintf(fid, "CNODE_{M} robustness: %f with average computation time: %f seconds.\n\n", sum(rMid.rob_ode)/rMid.numT, rMid.timeT/rMid.numT);
    fprintf(fid, "Control System - ACC\n");
    fprintf(fid, "ACC Computation time: %f seconds\n\n\n", rACC.rT);
    
    % Section 4.3
    fprintf(fid, "-------------------------------------\n");
    fprintf(fid, "SECTION 4.3 - RNNs \n");
    fprintf(fid, "-------------------------------------\n\n");
    fprintf(fid, "R_2_0 is robust in %i out of 20 experiments.\n", compute_robustness(rR2.rb1));
    fprintf(fid, "R_4_4 is robust in %i out of 20 experiments.\n", compute_robustness(rR4.rb1));
    fprintf(fid, "R_8_0 is robust in %i out of 20 experiments.\n\n\n", compute_robustness(rR8.rb1));
    
    % Section 4.4
    fprintf(fid, "-------------------------------------\n");
    fprintf(fid, "SECTION 4.4 - Semantic Segmentation \n");
    fprintf(fid, "-------------------------------------\n\n");
    fprintf(fid, "Dilated SSNN \n");
    fprintf(fid, "    Robustness: %f %%, Sensitivity: %f, IoU: %f, and computation time: %f seconds.\n", rDil.rv*100, rDil.rs, rDil.riou, rDil.t);
    fprintf(fid, "Transposed SSNN \n");
    fprintf(fid, "    Robustness: %f %%, Sensitivity: %f, IoU: %f, and computation time: %f seconds.\n", rTrans.rv*100, rTrans.rs, rTrans.riou, rTrans.t);
    fclose(fid);

end

% Create txt file with experiements from 4.2-4.4
function other_compare()

    % TABLE 3 - RL, tllverify, oval 21
    oval = load("NNV_vs_MATLAB/oval21/results_oval21.mat");
    ovalT = sum(oval.res(:,2))/length(oval.res(:,2));
    ovalS = sum(oval.res(:,1) == 0);
    ovalU = sum(oval.res(:,1) == 1);
    rl = load("NNV_vs_MATLAB/rl_benchmarks/results_rl.mat");
    rlTm = sum(rl.res(:,2))/length(rl.res(:,2));
    rlSm = sum(rl.res(:,1) == 0);
    rlUm = sum(rl.res(:,1) == 1);
    rlT = sum(rl.res(:,4))/length(rl.res(:,4));
    rlS = sum(rl.res(:,3) == 0);
    rlU = sum(rl.res(:,3) == 1);
    tll = load("NNV_vs_MATLAB/tllverify/results_tllverify.mat");
    tllTm = sum(tll.res(:,2))/length(tll.res(:,2));
    tllSm = sum(tll.res(:,1) == 0);
    tllUm = sum(tll.res(:,1) == 1);
    tllT = sum(tll.res(:,4))/length(tll.res(:,4));
    tllS = sum(tll.res(:,3) == 0);
    tllU = sum(tll.res(:,3) == 1);
    if is_codeocean()
        fid = fopen([path_results_codeocean, 'Table_3.txt'],'w'); 
    else
        fid = fopen('Table_3.txt','w'); 
    end
    fprintf(fid, '-----------------------------------------------------------\n');
    fprintf(fid, "  MATLAB  \n");
    fprintf(fid, '-----------------------------------------------------------\n');
    fprintf(fid, "RL Benchmarks (50) \n");
    fprintf(fid, "    SAT: %i, UNSAT: %i, time (s): %.4g\n", rlSm, rlUm, rlTm);
    fprintf(fid, "Tllverify Benchmarks (10) \n");
    fprintf(fid, "    SAT: %i, UNSAT: %i, time (s): %.4g\n", tllSm, tllUm, tllTm);
    fprintf(fid, "Oval21 Benchmarks (30) \n");
    fprintf(fid, "    Not supported \n");
    fprintf(fid, '-----------------------------------------------------------\n');
    fprintf(fid, "  NNV  \n");
    fprintf(fid, '-----------------------------------------------------------\n');
    fprintf(fid, "RL Benchmarks (50) \n");
    fprintf(fid, "    SAT: %i, UNSAT: %i, time (s): %.4g\n", rlS, rlU, rlT);
    fprintf(fid, "Tllverify Benchmarks (10) \n");
    fprintf(fid, "    SAT: %i, UNSAT: %i, time (s): %.4g\n", tllS, tllU, tllT);
    fprintf(fid, "Oval21 Benchmarks (30) \n");
    fprintf(fid, "    SAT: %i, UNSAT: %i, time (s): %.4g\n", ovalS, ovalU, ovalT);
end

% Complete RE package 
function pLong()
    % First do the short script (all - comparison)
    pShort();
    % Comparison to matlab
    % ACAS Xu (Table 2)
    acas_compare();
    % Other (Table 3)
    other_compare();
end
