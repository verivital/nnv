function run_cav23()

    % Reproduce all CAV23 experiments

    % Turn off figure display
    set(0,'DefaultFigureVisible','off');

    % Turn warnings off
    warning('off','all');
        
    %% 1) Neural ODEs
    disp("Running neural ODE experiments...");
    cd NeuralODEs;
    run_neuralode;
    cd ..;
    disp("Neural ODE experiments finished.");

    %% 2) NNV_vs_MATLAB

    disp("Running NNV vs MATLAB comparison experiments...");

    cd NNV_vs_MATLAB;

    disp("Running ACAS Xu experiments...");
    % ACAS Xu
    cd acas;
    verifyAll;
    cd ..;

    disp("Running OVAL21 experiments...");
    % OVAL 21
    cd oval21;
    verifyAll;
    cd ..;

    disp("Running RL Benchmarks experiments...");
    % RL benchmarks
    cd rl_benchmarks;
    verifyAll;
    cd ..;
    
    disp("Running TLLverify experiments...");
    % TLLverify
    cd tllverify;
    verifyAll;
    cd ..;
    
    cd ..;

    disp("NNV vs MATLAB comparison is finished.");

    %% 3) Segmentation

    disp("Running Semantic Segmentation examples...");

    cd Segmentation;
    example_dilated;
    example_transposed;
    cd ..;

    disp("Semantic Segmentation examples are finished.");

    %% 4) RNN

    disp("Running RNNs experiments...");

    cd RNN;
    run_rnn_all();
    create_figure();
    cd ..

    disp("RNN examples are finished.");

    %% Print results
    
    disp("Generating manuscript figures and tables...");
    paper_results('all');

    disp("SUCCESFULLY FINISHED ALL EXPERIMENTS.")

end

