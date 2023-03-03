function run_rnn_all()

    % Run all verification experiments for RNNs
    cd N_2_0;
    verify_N2_0;
    cd ..;

    cd N_4_4;
    verify_N4_4;
    cd ..;

    cd N_8_0
    verify_N8_0;
    cd ..;
end

