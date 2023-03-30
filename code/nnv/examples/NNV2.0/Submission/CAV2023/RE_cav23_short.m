function RE_cav23_short()

    % Reproduce a subset of CAV23 experiments

    % Turn off figure display
    set(0,'DefaultFigureVisible','off');

    % Turn warnings off
    warning('off','all');
        
    %% 1) Neural ODEs

    cd NeuralODEs;
    run_neuralode;
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
    
    paper_results('short');

end

