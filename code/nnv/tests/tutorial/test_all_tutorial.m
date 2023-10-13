% Adding all the tutorial scripts to the tests to ensure it always works

cd ../../examples/Tutorial/NN/'ACAS Xu'/;

%% 1) NN: ACAS Xu (verify_onnx_vnnlib)
verify_onnx_vnnlib;

%% 2) NN: GTSRB (verify_falsify_1)
cd ../GTSRB;
verify_falsify_1;

%% 3) NN: GTSRB (verify_robust_13)
verify_robust_13;

%% 4) NN: GTSRB (verify_robust_27)
verify_robust_27;

%% 5) NN: GTSRB (verify_unknown_1)
verify_unknown_1;

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

%% 12) NNCS: AEBS (reach_full_brake)
reach_full_brake;

%% 13) NNCS: Inverted Pendulum
cd ../InvertedPendulum;
reach_invP;

%% 14) Other: load models
cd ../../other;
load_models;

%% 15) Other: set representations
set_representations;
