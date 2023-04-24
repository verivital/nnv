function run_rnn_all()

    % Run all verification experiments for RNNs
    disp("Running RNN N_2_0 analysis...");
    cd N_2_0;
    verify_N2_0;
    cd ..;

    disp("Running RNN N_4_4 analysis...");
    cd N_4_4;
    verify_N4_4;
    cd ..;

    disp("Running RNN N_8_0 analysis...")
    cd N_8_0
    verify_N8_0;
    cd ..;
    
end

