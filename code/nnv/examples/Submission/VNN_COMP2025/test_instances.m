vnncomp_path = "C:\Users\diego\Documents\Research\vnncomp2025_benchmarks\benchmarks\";

% Things look better than last year, let's make sure we have no penalties this time
% Can we support any other benchmarks? what are the errors we are getting
% in some of them?


%% acasxu

disp("Running acas xu...")

acas_path = vnncomp_path + "acasxu_2023/";

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
    run_vnncomp_instance("nn4sys",onnx,vnnlib,"nn4sys_results_" + string(i)+".txt");
end

% need to update the support for these vnnlib properties, should be pretty quick

%% dist_shift

disp("Running dist_shift...")

dist_path = vnncomp_path + "dist_shift_2023/";

dist_instances = ["onnx/mnist_concat.onnx" ,"vnnlib/index6731_delta0.13.vnnlib";...
    "onnx/mnist_concat.onnx" ,"vnnlib/index7_delta0.13.vnnlib";...
    "onnx/mnist_concat.onnx" ,"vnnlib/index9830_delta0.13.vnnlib";...
    ]; % something is failing with these ones, let's take a look at that

% Run verification for dist_shift
for i=1:length(dist_instances)
    onnx = dist_path + dist_instances(i,1);
    vnnlib = dist_path + dist_instances(i,2);
    run_vnncomp_instance("dist_shift",onnx,vnnlib,"dist_results_" + string(i)+".txt");
end


%% collins_rul

disp("Running collins_rul..")

rul_path = vnncomp_path + "collins_rul_cnn_2022/";

rul_instances = ["onnx/NN_rul_small_window_20.onnx" ,"vnnlib/robustness_2perturbations_delta5_epsilon10_w20.vnnlib";...
    "onnx/NN_rul_small_window_20.onnx" ,"vnnlib/monotonicity_CI_shift5_w20.vnnlib";...
    "onnx/NN_rul_small_window_20.onnx" ,"vnnlib/if_then_5levels_w20.vnnlib";...
    "onnx/NN_rul_full_window_20.onnx", "vnnlib/robustness_2perturbations_delta5_epsilon10_w20.vnnlib";...
    "onnx/NN_rul_full_window_40.onnx" ,"vnnlib/robustness_2perturbations_delta5_epsilon10_w40.vnnlib"];

% Run verification for collins_rul
for i=1:length(rul_instances)
    onnx = rul_path + rul_instances(i,1);
    vnnlib = rul_path + rul_instances(i,2);
    run_vnncomp_instance("collins_rul",onnx,vnnlib,"collins_rul_results_" + string(i)+".txt");
end

%% cgan

disp("Running cgan..")

cgan_path = vnncomp_path + "cgan_2023/";

cgan_instances = ["onnx/cGAN_imgSz32_nCh_1.onnx", "vnnlib/cGAN_imgSz32_nCh_1_prop_0_input_eps_0.010_output_eps_0.015.vnnlib";...
    "onnx/cGAN_imgSz32_nCh_3.onnx", "vnnlib/cGAN_imgSz32_nCh_3_prop_0_input_eps_0.015_output_eps_0.020.vnnlib";...
    "onnx/cGAN_imgSz64_nCh_1.onnx", "vnnlib/cGAN_imgSz64_nCh_1_prop_1_input_eps_0.005_output_eps_0.010.vnnlib";...
    "onnx/cGAN_imgSz64_nCh_3.onnx", "vnnlib/cGAN_imgSz64_nCh_3_prop_2_input_eps_0.005_output_eps_0.010.vnnlib";...
    % "onnx/cGAN_imgSz32_nCh_3_nonlinear_activations.onnx", "vnnlib/cGAN_imgSz32_nCh_3_nonlinear_activations_prop_0_input_eps_0.010_output_eps_0.015.vnnlib";... % conv2D error (data and filter have different data types)
    "onnx/cGAN_imgSz32_nCh_1_transposedConvPadding_1.onnx", "vnnlib/cGAN_imgSz32_nCh_1_transposedConvPadding_1_prop_0_input_eps_0.010_output_eps_0.015.vnnlib";...
    % "onnx/cGAN_imgSz32_nCh_3_upsample.onnx", "vnnlib/cGAN_imgSz32_nCh_3_upsample_prop_0_input_eps_0.010_output_eps_0.015.vnnlib";... % reshape error reachability
    % "onnx/cGAN_imgSz32_nCh_3_small_transformer.onnx", "vnnlib/cGAN_imgSz32_nCh_3_small_transformer_prop_0_input_eps_0.005_output_eps_0.010.vnnlib";... % onnx loading error
];

% Run verification for cgan
for i=1:length(cgan_instances)
    onnx = cgan_path + cgan_instances(i,1);
    vnnlib = cgan_path + cgan_instances(i,2);
    run_vnncomp_instance("cgan",onnx,vnnlib,"cgan_results_" + string(i)+".txt");
end


%% tllverify

disp("Running tllverify..")

tll_path = vnncomp_path + "tllverifybench_2023/";

tll_instances = ["onnx/tllBench_n=2_N=M=8_m=1_instance_0_0.onnx","vnnlib/property_N=8_0.vnnlib";...
    "onnx/tllBench_n=2_N=M=16_m=1_instance_1_0.onnx","vnnlib/property_N=16_0.vnnlib";...
    "onnx/tllBench_n=2_N=M=24_m=1_instance_2_2.onnx","vnnlib/property_N=24_2.vnnlib";...
    "onnx/tllBench_n=2_N=M=32_m=1_instance_3_0.onnx","vnnlib/property_N=32_0.vnnlib";...
    "onnx/tllBench_n=2_N=M=48_m=1_instance_5_3.onnx","vnnlib/property_N=48_3.vnnlib";...
    "onnx/tllBench_n=2_N=M=56_m=1_instance_6_0.onnx","vnnlib/property_N=56_0.vnnlib";...
    % "onnx/tllBench_n=2_N=M=64_m=1_instance_7_0.onnx","vnnlib/property_N=64_0.vnnlib";...
    ]; % last one not finishing after a loooong time

% Run verification for tllverify
for i=1:length(tll_instances)
    onnx = tll_path + tll_instances(i,1);
    vnnlib = tll_path + tll_instances(i,2);
    run_vnncomp2024_instance("tllverifybench",onnx,vnnlib,"tllverify_results_" + string(i)+".txt");
end

%% metaroom

disp("Running metaroom..")

metaroom_path = vnncomp_path + "metaroom_2023/";

% We can probably do exact analysis from the beginning on all instances here
metaroom_instances = ["onnx/6cnn_ry_20_3_no_custom_OP.onnx", "vnnlib/spec_idx_109_eps_0.00000436.vnnlib";...
    "onnx/4cnn_ry_117_19_no_custom_OP.onnx", "vnnlib/spec_idx_3_eps_0.00000436.vnnlib";...
    "onnx/6cnn_ry_65_10_no_custom_OP.onnx", "vnnlib/spec_idx_133_eps_0.00000436.vnnlib";...
    "onnx/6cnn_ry_104_17_no_custom_OP.onnx", "vnnlib/spec_idx_96_eps_0.00000436.vnnlib";...
    "onnx/6cnn_ry_94_15_no_custom_OP.onnx", "vnnlib/spec_idx_149_eps_0.00000436.vnnlib";...
    "onnx/4cnn_ry_19_3_no_custom_OP.onnx", "vnnlib/spec_idx_5_eps_0.00000436.vnnlib";...
    ];

% Run verification for metaroom
for i=1:length(metaroom_instances)
    onnx = metaroom_path + metaroom_instances(i,1);
    vnnlib = metaroom_path + metaroom_instances(i,2);
    run_vnncomp_instance("metaroom",onnx,vnnlib,"metaroom_results_" + string(i)+".txt");
end


%% safeNLP

disp("Running safeNLP..")

safeNLP_path = vnncomp_path + "safenlp_2024/";

% All are unsat in about a second, we can prob do exact on all as well

% We can probably do exact analysis from the beginning on all instances here
safeNLP_instances = ["onnx/ruarobot.onnx", "vnnlib/ruarobot/prop_8_perturbations.vnnlib";...
    "onnx/ruarobot.onnx", "vnnlib/ruarobot/prop_1239_perturbations.vnnlib";...
    "onnx/ruarobot.onnx", "vnnlib/ruarobot/prop_6012_perturbations.vnnlib";...
    "onnx/medical.onnx", "vnnlib/medical/prop_0_perturbations.vnnlib";...
    "onnx/medical.onnx", "vnnlib/medical/prop_5_perturbations.vnnlib";...
    "onnx/medical.onnx", "vnnlib/medical/prop_53_perturbations.vnnlib";...
    ];

% Run verification for safeNLP
for i=1:length(safeNLP_instances)
    onnx = safeNLP_path + safeNLP_instances(i,1);
    vnnlib = safeNLP_path + safeNLP_instances(i,2);
    run_vnncomp_instance("safeNLP",onnx,vnnlib,"safeNLP_results_" + string(i)+".txt");
end


%% Cora

disp("Running Cora benchmark..")

cora_path = vnncomp_path + "cora_2024/";

% Some sat, most unsat, fairly fast using exact, can go with that as well
% Some unknowns... How is this possible with exact analysis?

% We can probably do exact analysis from the beginning on all instances here
cora_instances = ["nns/mnist-point.onnx", "benchmark-files/mnist-img0.vnnlib";...
    "nns/mnist-trades.onnx", "benchmark-files/mnist-img1.vnnlib";...
    "nns/mnist-set.onnx", "benchmark-files/mnist-img20.vnnlib";...
    "nns/svhn-point.onnx", "benchmark-files/svhn-img140.vnnlib";...
    "nns/svhn-set.onnx", "benchmark-files/svhn-img323.vnnlib";...
    "nns/svhn-trades.onnx", "benchmark-files/svhn-img472.vnnlib";...
    "nns/cifar10-set.onnx", "benchmark-files/cifar10-img347.vnnlib";...
    "nns/cifar10-point.onnx", "benchmark-files/cifar10-img353.vnnlib";...
    "nns/cifar10-trades.onnx", "benchmark-files/cifar10-img266.vnnlib"...
    ];

% Run verification for CORA
for i=1:length(cora_instances)
    onnx = cora_path + cora_instances(i,1);
    vnnlib = cora_path + cora_instances(i,2);
    run_vnncomp_instance("cora",onnx,vnnlib,"cora_results_" + string(i)+".txt");
end


%% Cifar100

disp("Running Cifar-100..")

cifar100_path = vnncomp_path + "cifar100_2024/";

% no vnnlib properties provided in this year's repo
cifar100_instances = ["onnx/CIFAR100_resnet_medium.onnx", ""; ...
    "onnx/CIFAR100_resnet_large.onnx", "";...
    ];

% Run verification for cifar100
for i=1:length(cifar100_instances)
    onnx = cifar100_path + cifar100_instances(i,1);
    vnnlib = cifar100_path + cifar100_instances(i,2);
    run_vnncomp_instance("cifar100",onnx,vnnlib,"cifar100_results_" + string(i)+".txt");
end


%% tinyimagenet

disp("Running tinyimagenet..")

tinyimagenet_path = vnncomp_path + "tinyimagenet_2024/";

% no vnnlib properties provided in this year's repo
% (Some may be available in 2022 repo)
tinyimagenet_instances = ["onnx/TinyImageNet_resnet_medium.onnx", ""; ...
    ];

% Run verification for tinyimagenet
for i=1:length(tinyimagenet_instances)
    onnx = tinyimagenet_path + tinyimagenet_instances(i,1);
    vnnlib = tinyimagenet_path + tinyimagenet_instances(i,2);
    run_vnncomp_instance("tinyimagenet",onnx,vnnlib,"tinyimagenet_results_" + string(i)+".txt");
end

%% TODO
% Add all other benchmarks
