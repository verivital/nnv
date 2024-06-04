%% Run some verification instances from the competition

% path to benchmarks
vnncomp_path = "/home/manzand/Documents/MATLAB/vnncomp2023_benchmarks/benchmarks/"; % change to whatever local path vnncomp benchmarks are stored
% vnncomp_path = "/home/dieman95/Documents/MATLAB/vnncomp2023_benchmarks/benchmarks/";


%% acasxu

disp("Running acas xu...")

acas_path = vnncomp_path + "acasxu/";

acas_instances = [...
    "onnx/ACASXU_run2a_1_1_batch_2000.onnx" , "vnnlib/prop_1.vnnlib";...
    "onnx/ACASXU_run2a_2_3_batch_2000.onnx","vnnlib/prop_2.vnnlib";...
    "onnx/ACASXU_run2a_3_4_batch_2000.onnx","vnnlib/prop_3.vnnlib";...
    "onnx/ACASXU_run2a_2_5_batch_2000.onnx","vnnlib/prop_4.vnnlib";...
    "onnx/ACASXU_run2a_1_1_batch_2000.onnx","vnnlib/prop_5.vnnlib";...
    "onnx/ACASXU_run2a_1_1_batch_2000.onnx","vnnlib/prop_6.vnnlib";...
    "onnx/ACASXU_run2a_1_9_batch_2000.onnx","vnnlib/prop_7.vnnlib";...
    "onnx/ACASXU_run2a_2_9_batch_2000.onnx","vnnlib/prop_8.vnnlib";...
    "onnx/ACASXU_run2a_3_3_batch_2000.onnx","vnnlib/prop_9.vnnlib";...
    "onnx/ACASXU_run2a_4_5_batch_2000.onnx","vnnlib/prop_10.vnnlib";...
    ];

% Run verification for acas 
for i=1:length(acas_instances)
    disp(i);
    onnx = acas_path + acas_instances(i,1);
    vnnlib = acas_path + acas_instances(i,2);
    run_vnncomp2024_instance("acasxu",onnx,vnnlib,"acas_results_" + string(i)+".txt");
end

%% test 

disp("Running test benchmark...")

test_path = vnncomp_path + "test/";

test_instances = [...
    "test_sat.onnx" , "test_prop.vnnlib";...
    "test_unsat.onnx" , "test_prop.vnnlib"];

% Run verification for test
for i=1:length(test_instances)
    onnx = test_path + test_instances(i,1);
    vnnlib = test_path + test_instances(i,2);
    run_vnncomp2024_instance("test",onnx,vnnlib,"test_results_" + string(i)+".txt");
end


%% nn4sys

disp("Running nn4sys...")

nn4sys_path = vnncomp_path + "nn4sys/";

nn4sys_instances = [... % all other networks are not supported...
    "onnx/lindex.onnx","vnnlib/lindex_1.vnnlib";...
    "onnx/lindex_deep.onnx", "vnnlib/lindex_200.vnnlib";...
    ];

% Run verification for nn4sys 
for i=1:length(nn4sys_instances)
    onnx = nn4sys_path + nn4sys_instances(i,1);
    vnnlib = nn4sys_path + nn4sys_instances(i,2);
    run_vnncomp2024_instance("nn4sys",onnx,vnnlib,"nn4sys_results_" + string(i)+".txt");
end

% need to update the support for these vnnlib properties, should be pretty quick

%% dist_shift

disp("Running dist_shift...")

dist_path = vnncomp_path + "dist_shift/";

dist_instances = ["onnx/mnist_concat.onnx" ,"vnnlib/index7165_delta0.13.vnnlib"];

% Run verification for nn4sys
onnx = dist_path + dist_instances(1,1);
vnnlib = dist_path + dist_instances(1,2);
run_vnncomp2024_instance("dist_shift",onnx,vnnlib,"dist_results_" + string(i)+".txt");


%% collins_rul

disp("Running collins_rul..")

rul_path = vnncomp_path + "collins_rul_cnn/";

rul_instances = ["onnx/NN_rul_small_window_20.onnx" ,"vnnlib/robustness_2perturbations_delta5_epsilon10_w20.vnnlib";...
    "onnx/NN_rul_small_window_20.onnx" ,"vnnlib/monotonicity_CI_shift5_w20.vnnlib";...
    "onnx/NN_rul_small_window_20.onnx" ,"vnnlib/if_then_5levels_w20.vnnlib";...
    "onnx/NN_rul_full_window_20.onnx", "vnnlib/robustness_2perturbations_delta5_epsilon10_w20.vnnlib";...
    "onnx/NN_rul_full_window_40.onnx" ,"vnnlib/robustness_2perturbations_delta5_epsilon10_w40.vnnlib"];

% Run verification for collins_rul
for i=1:length(rul_instances)
    onnx = rul_path + rul_instances(i,1);
    vnnlib = rul_path + rul_instances(i,2);
    run_vnncomp2024_instance("collins_rul",onnx,vnnlib,"collins_rul_results_" + string(i)+".txt");
end

%% cgan

disp("Running cgan..")

cgan_path = vnncomp_path + "cgan/";

cgan_instances = ["onnx/cGAN_imgSz64_nCh_3.onnx" ,"vnnlib/cGAN_imgSz64_nCh_3_prop_3_input_eps_0.010_output_eps_0.015.vnnlib";...
    "onnx/cGAN_imgSz32_nCh_3_nonlinear_activations.onnx","vnnlib/cGAN_imgSz32_nCh_3_nonlinear_activations_prop_0_input_eps_0.015_output_eps_0.020.vnnlib";...
    "onnx/cGAN_imgSz32_nCh_1_transposedConvPadding_1.onnx" ,"vnnlib/cGAN_imgSz32_nCh_1_transposedConvPadding_1_prop_0_input_eps_0.015_output_eps_0.020.vnnlib";...
    "onnx/cGAN_imgSz32_nCh_3_upsample.onnx","vnnlib/cGAN_imgSz32_nCh_3_upsample_prop_0_input_eps_0.015_output_eps_0.020.vnnlib"];

% Run verification for nn4sys
for i=1:length(cgan_instances)
    onnx = cgan_path + cgan_instances(i,1);
    vnnlib = cgan_path + cgan_instances(i,2);
    run_vnncomp2024_instance("cgan",onnx,vnnlib,"cgan_results_" + string(i)+".txt");
end


%% ml4acopf

disp("Running ml4acopf..")

ml4_path = vnncomp_path + "ml4acopf/";

ml4_instances = ["onnx/118_ieee_ml4acopf.onnx","vnnlib/118_ieee_prop2.vnnlib";...
    "onnx/14_ieee_ml4acopf.onnx","vnnlib/14_ieee_prop9.vnnlib";...
    "onnx/300_ieee_ml4acopf.onnx","vnnlib/300_ieee_prop3.vnnlib"];

% Run verification for nn4sys
for i=1:length(ml4_instances)
    onnx = ml4_path + ml4_instances(i,1);
    vnnlib = ml4_path + ml4_instances(i,2);
    run_vnncomp2024_instance("ml4acopf",onnx,vnnlib,"ml4_results_" + string(i)+".txt");
end


%% tllverify

disp("Running tllverify..")

tll_path = vnncomp_path + "tllverifybench/";

tll_instances = ["onnx/tllBench_n=2_N=M=8_m=1_instance_0_0.onnx","vnnlib/property_N=8_0.vnnlib";...
    "onnx/tllBench_n=2_N=M=16_m=1_instance_1_0.onnx","vnnlib/property_N=16_0.vnnlib";...
    "onnx/tllBench_n=2_N=M=24_m=1_instance_2_2.onnx","vnnlib/property_N=24_2.vnnlib";...
    "onnx/tllBench_n=2_N=M=32_m=1_instance_3_0.onnx","vnnlib/property_N=32_0.vnnlib";...
    "onnx/tllBench_n=2_N=M=48_m=1_instance_5_3.onnx","vnnlib/property_N=48_3.vnnlib";...
    "onnx/tllBench_n=2_N=M=56_m=1_instance_6_0.onnx","vnnlib/property_N=56_0.vnnlib";...
    "onnx/tllBench_n=2_N=M=64_m=1_instance_7_0.onnx","vnnlib/property_N=64_0.vnnlib"];

% Run verification for nn4sys
for i=1:length(tll_instances)
    onnx = tll_path + tll_instances(i,1);
    vnnlib = tll_path + tll_instances(i,2);
    run_vnncomp2024_instance("tllverifybench",onnx,vnnlib,"tllverify_results_" + string(i)+".txt");
end


%% TODO
% Add benchmarks from 2024
% update nn4sys and dist_shift as they have been updated