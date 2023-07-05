%% Run some verification instances from the competition

% path to benchmarks
vnncomp_path = "/home/manzand/Documents/MATLAB/vnncomp2023_benchmarks/benchmarks/";
% vnncomp_path = "/home/dieman95/Documents/MATLAB/vnncomp2023_benchmarks/benchmarks/";

% Go through some of the instances for every benchmark

%% acasxu
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
    onnx = acas_path + acas_instances(i,1);
    vnnlib = acas_path + acas_instances(i,2);
    run_vnncomp_instance("acasxu",onnx,vnnlib,"acas_results_" + string(i)+".txt");
end

%% ViT

vit_path = vnncomp_path + "vit/";

vit_instances = [...
    "onnx/ibp_3_3_8.onnx" , "vnnlib/ibp_3_3_8_3850.vnnlib";...
    "onnx/pgd_2_3_16.onnx","vnnlib/pgd_2_3_16_1744.vnnlib";...
    ];

% Run verification for vit
for i=1:length(vit_instances)
    onnx = vit_path + vit_instances(i,1);
    vnnlib = vit_path + vit_instances(i,2);
    run_vnncomp_instance("vit",onnx,vnnlib,"vit_results_" + string(i)+".txt");
end


%% nn4sys

nn4sys_path = vnncomp_path + "nn4sys/";

nn4sys_instances = [... % no instances found with pensieve_mid_parallel
    "onnx/pensieve_big_simple.onnx", "vnnlib/pensieve_simple_0.vnnlib";...
    "onnx/pensieve_small_parallel.onnx","vnnlib/pensieve_parallel_1.vnnlib";...
    "onnx/pensieve_mid_simple.onnx" ," vnnlib/pensieve_simple_2.vnnlib";...
    "onnx/pensieve_big_parallel.onnx", "vnnlib/pensieve_parallel_35.vnnlib";...
    "onnx/pensieve_small_simple.onnx", "vnnlib/pensieve_simple_24.vnnlib";...
    "onnx/lindex.onnx","vnnlib/lindex_1.vnnlib";...
    "onnx/lindex_deep.onnx", "vnnlib/lindex_9000.vnnlib";...
    "onnx/mscn_128d.onnx" , "vnnlib/cardinality_0_1_128.vnnlib";...
    "onnx/mscn_128d_dual.onnx" , "vnnlib/cardinality_1_240_128_dual.vnnlib";...
    "onnx/mscn_2048d.onnx" ,"vnnlib/cardinality_0_500_2048.vnnlib";...
    "onnx/mscn_2048d_dual.onnx" ,"vnnlib/cardinality_1_1_2048_dual.vnnlib";...
    ];

% Run verification for nn4sys
for i=1:length(nn4sys_instances)
    onnx = nn4sys_path + nn4sys_instances(i,1);
    vnnlib = nn4sys_path + nn4sys_instances(i,2);
    run_vnncomp_instance("nn4sys",onnx,vnnlib,"nn4sys_results_" + string(i)+".txt");
end


%% dist_shift

dist_path = vnncomp_path + "dist_shift/";

dist_instances = ["onnx/mnist_concat.onnx" ,"vnnlib/index7165_delta0.13.vnnlib"];

% Run verification for nn4sys
for i=1:length(dist_instances)
    onnx = dist_path + dist_instances(i,1);
    vnnlib = dist_path + dist_instances(i,2);
    run_vnncomp_instance("dist_shift",onnx,vnnlib,"dist_results_" + string(i)+".txt");
end


%% traffic_sign
% we should be able to support it, but matlab does not, so we would have to create the models manually...
% can load most info with Keras importer

%% collins_rul

rul_path = vnncomp_path + "collins_rul_cnn/";

rul_instances = ["onnx/NN_rul_small_window_20.onnx" ,"vnnlib/robustness_2perturbations_delta5_epsilon10_w20.vnnlib";...
    "onnx/NN_rul_small_window_20.onnx" ,"vnnlib/monotonicity_CI_shift5_w20.vnnlib";...
    "onnx/NN_rul_small_window_20.onnx" ,"vnnlib/if_then_5levels_w20.vnnlib";...
    "onnx/NN_rul_full_window_20.onnx", "vnnlib/robustness_2perturbations_delta5_epsilon10_w20.vnnlib";...
    "onnx/NN_rul_full_window_40.onnx" ,"vnnlib/robustness_2perturbations_delta5_epsilon10_w40.vnnlib"];

% Run verification for nn4sys
for i=1:length(rul_instances)
    onnx = rul_path + rul_instances(i,1);
    vnnlib = rul_path + rul_instances(i,2);
    run_vnncomp_instance("collins_rul",onnx,vnnlib,"collins_rul_results_" + string(i)+".txt");
end

%% cgan

cgan_path = vnncomp_path + "cgan/";

cgan_instances = ["onnx/cGAN_imgSz64_nCh_3.onnx" ,"vnnlib/cGAN_imgSz64_nCh_3_prop_3_input_eps_0.010_output_eps_0.015.vnnlib";...
    "onnx/cGAN_imgSz32_nCh_3_nonlinear_activations.onnx","vnnlib/cGAN_imgSz32_nCh_3_nonlinear_activations_prop_0_input_eps_0.015_output_eps_0.020.vnnlib";...
    "onnx/cGAN_imgSz32_nCh_1_transposedConvPadding_1.onnx" ,"vnnlib/cGAN_imgSz32_nCh_1_transposedConvPadding_1_prop_0_input_eps_0.015_output_eps_0.020.vnnlib";...
    "onnx/cGAN_imgSz32_nCh_3_upsample.onnx","vnnlib/cGAN_imgSz32_nCh_3_upsample_prop_0_input_eps_0.015_output_eps_0.020.vnnlib";...
    "onnx/cGAN_imgSz32_nCh_3_small_transformer.onnx","vnnlib/cGAN_imgSz32_nCh_3_small_transformer_prop_1_input_eps_0.010_output_eps_0.015.vnnlib"];

% Run verification for nn4sys
for i=1:length(cgan_instances)
    onnx = cgan_path + cgan_instances(i,1);
    vnnlib = cgan_path + cgan_instances(i,2);
    run_vnncomp_instance("cgan",onnx,vnnlib,"cgan_results_" + string(i)+".txt");
end

%% vggnet16

vgg_path = vnncomp_path + "vggnet16/";

vgg_instances = ["onnx/vgg16-7.onnx","vnnlib/spec0_screw.vnnlib"];

% Run verification for nn4sys
for i=1:length(vgg_instances)
    onnx = vgg_path + vgg_instances(i,1);
    vnnlib = vgg_path + vgg_instances(i,2);
    run_vnncomp_instance("vggnet16",onnx,vnnlib,"vgg_results_" + string(i)+".txt");
end

%% ml4acopf

ml4_path = vnncomp_path + "ml4acopf/";

ml4_instances = ["onnx/118_ieee_ml4acopf.onnx","vnnlib/118_ieee_prop2.vnnlib";...
    "onnx/14_ieee_ml4acopf.onnx","vnnlib/14_ieee_prop9.vnnlib";...
    "onnx/300_ieee_ml4acopf.onnx","vnnlib/300_ieee_prop3.vnnlib"];

% Run verification for nn4sys
for i=1:length(ml4_instances)
    onnx = ml4_path + ml4_instances(i,1);
    vnnlib = ml4_path + ml4_instances(i,2);
    run_vnncomp_instance("ml4acopf",onnx,vnnlib,"ml4_results_" + string(i)+".txt");
end


%% tllverify

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
    run_vnncomp_instance("tllverifybench",onnx,vnnlib,"tllverify_results_" + string(i)+".txt");
end

