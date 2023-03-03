function mnist_table()

    %% Create MNIST results table

    % Run 1: inf 0.5
    cnm = load("cnn_medium_nnv_inf_0.5.mat");
    cns = load("cnn_tiny_nnv_inf_0.5.mat");

    % Create tables
    robs =  [cns.rob/cns.numT;   cnm.rob/cnm.numT];
    timeA = [cns.timeT/cns.numT; cnm.timeT/cnm.numT;];
    
    % names = ["CNODEsmall";"CNODEmid"];
    % colmns = {'Robust' ;'Time (s)' 
    T = table(robs,timeA);
    
    % T.Properties.VariableNames = colmns;
    if is_codeocean
        table2latex(T,'/results/logs/table_minist.tex')
    else
        table2latex(T,'table_minist.tex')
    end

end