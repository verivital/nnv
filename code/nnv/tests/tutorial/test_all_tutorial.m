% Adding all the tutorial scripts to the tests to ensure it always works

cd ../../examples/Tutorial/NN/'ACAS Xu'/;

%% 1) NN: ACAS Xu (verify_onnx_vnnlib)
if ~is_github_actions() % importers (support packages) are not installed
    verify_onnx_vnnlib;
end


%% 3) NN: GTSRB (verify_robust_13)
cd ../GTSRB;
verify_robust_13;

%% 6) NN: MNIST (input set examples)
cd ../MNIST;
input_set_examples;

%% 7) NN: MNIST (verify)
verify;

%% 8) NN: MNIST (verify_fc)
verify_fc;

%% 9) NN: Segmentation
cd ../Segmentation;
verify_m2nist;

%% 10) NNCS: ACC
cd ../../NNCS/ACC/Verification;
verify;

%% 11) NNCS: AEBS (reach)
cd ../../AEBS;
reach;

%% 13) NNCS: Inverted Pendulum
cd ../InvertedPendulum;
reach_invP;

%% 14) Other: load models
cd ../../other;
if ~is_github_actions() % importers (support packages) are not installed
    load_models;
end

%% 15) Other: set representations
set_representations;
