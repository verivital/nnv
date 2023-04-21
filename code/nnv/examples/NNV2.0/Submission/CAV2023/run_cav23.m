function run_cav23()

    % Reproduce all CAV23 experiments

    % Turn off figure display
    set(0,'DefaultFigureVisible','off');

    % Turn warnings off
    warning('off','all');
        
    %% 1) Neural ODEs

    cd NeuralODEs;
    run_neuralode;
    cd ..;

    %% 2) NNV_vs_MATLAB

    cd NNV_vs_MATLAB;

    % ACAS Xu
    cd acas;
    verifyAll;
    cd ..;

    % OVAL 21
    cd oval21;
    verifyAll;
    cd ..;

    % RL benchmarks
    cd rl_benchmarks;
    verifyAll;
    cd ..;
    
    % TLLverify
    cd tllverify;
    verifyAll;
    cd ..;
    
    cd ..;

    %% 3) Segmentation

    cd Segmentation;
    example_dilated;
    example_transposed;
    cd ..;

    %% 4) RNN

    cd RNN;
    run_rnn_all();
    create_figure();
    cd ..

    %% Print results
    
    paper_results('all');

end

