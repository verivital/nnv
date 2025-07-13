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
    try
        run_vnncomp_instance("acasxu",onnx,vnnlib,"acas_results_" + string(i)+".txt");
    catch ME
        warning("Failed")
        disp(onnx+"___"+vnnlib)
        warning(ME.message)
    end

end


% No errors


%% cctsdb_yolo

disp("Running cctsdb...")

cctsdb_path = vnncomp_path + "cctsdb_yolo_2023/";

cctsdb_instances = [...
    "onnx/patch-1.onnx" , "vnnlib/spec_onnx_patch-1_idx_00559_0.vnnlib";...
    "onnx/patch-3.onnx", "vnnlib/spec_onnx_patch-3_idx_00303_0.vnnlib";...
    ];

% Run verification for acas 
for i=1:length(cctsdb_instances)
    onnx = cctsdb_path + cctsdb_instances(i,1);
    vnnlib = cctsdb_path + cctsdb_instances(i,2);
    try
        run_vnncomp_instance("cctsdb_yolo",onnx,vnnlib,"cctsdb_results_" + string(i)+".txt");
    catch ME
        warning("Failed")
        disp(onnx+"___"+vnnlib)
        warning(ME.message)
    end
end

% Errors on both? None of them finished...
% Let's test this again

%% cersyve

disp("Running cersyve...")

cersyve_path = vnncomp_path + "cersyve/";

cersyve_instances = [...
    "onnx/unicycle_pretrain_con.onnx" , "vnnlib/prop_unicycle.vnnlib";...
    "onnx/robot_arm_finetune_inv.onnx", "vnnlib/prop_robot_arm.vnnlib";...
    ];

% Run verification for acas 
for i=1:length(cersyve_instances)
    onnx = cersyve_path + cersyve_instances(i,1);
    vnnlib = cersyve_path + cersyve_instances(i,2);
    try
        run_vnncomp_instance("cersyve",onnx,vnnlib,"cersyve_results_" + string(i)+".txt");
    catch ME
        warning("Failed")
        disp(onnx+"___"+vnnlib)
        warning(ME.message)
    end
end

% Both finished, no errors, but none finished to verify even on submission
% site


%% cgan

disp("Running cgan..")

cgan_path = vnncomp_path + "cgan_2023/";

cgan_instances = ["onnx/cGAN_imgSz32_nCh_1.onnx", "vnnlib/cGAN_imgSz32_nCh_1_prop_0_input_eps_0.010_output_eps_0.015.vnnlib";...
    "onnx/cGAN_imgSz32_nCh_3.onnx", "vnnlib/cGAN_imgSz32_nCh_3_prop_0_input_eps_0.010_output_eps_0.015.vnnlib";...
    "onnx/cGAN_imgSz32_nCh_3_nonlinear_activations.onnx", "vnnlib/cGAN_imgSz32_nCh_3_nonlinear_activations_prop_0_input_eps_0.015_output_eps_0.020.vnnlib";... 
    "onnx/cGAN_imgSz32_nCh_1_transposedConvPadding_1.onnx", "vnnlib/cGAN_imgSz32_nCh_1_transposedConvPadding_1_prop_0_input_eps_0.015_output_eps_0.020.vnnlib";...
    "onnx/cGAN_imgSz32_nCh_3_upsample.onnx", "vnnlib/cGAN_imgSz32_nCh_3_upsample_prop_0_input_eps_0.015_output_eps_0.020.vnnlib";... 
    "onnx/cGAN_imgSz32_nCh_3_small_transformer.onnx", "vnnlib/cGAN_imgSz32_nCh_3_small_transformer_prop_0_input_eps_0.010_output_eps_0.015.vnnlib";... 
];

% Run verification for cgan
for i=1:length(cgan_instances)
    onnx = cgan_path + cgan_instances(i,1);
    vnnlib = cgan_path + cgan_instances(i,2);
    try
        run_vnncomp_instance("cgan",onnx,vnnlib,"cgan_results_" + string(i)+".txt");
    catch ME
        warning("Failed")
        disp(onnx+"___"+vnnlib)
        warning(ME.message)
    end
end

% All good apparently, can even verify some as seen on submission site

%% Cifar100

disp("Running Cifar-100..")

cifar100_path = vnncomp_path + "cifar100_2024/";

% no vnnlib properties provided in this year's repo
cifar100_instances = ["onnx/CIFAR100_resnet_medium.onnx", "vnnlib/CIFAR100_resnet_medium_prop_idx_9697_sidx_8423_eps_0.0039.vnnlib"; ...
    "onnx/CIFAR100_resnet_large.onnx", "vnnlib/CIFAR100_resnet_large_prop_idx_2247_sidx_8763_eps_0.0039.vnnlib";...
    ];

% Run verification for cifar100
for i=1:length(cifar100_instances)
    onnx = cifar100_path + cifar100_instances(i,1);
    vnnlib = cifar100_path + cifar100_instances(i,2);
    try
        run_vnncomp_instance("cifar100",onnx,vnnlib,"cifar100_results_" + string(i)+".txt");
    catch ME
        warning("Failed")
        disp(onnx+"___"+vnnlib)
        warning(ME.message)
    end
end


% These will need to run on CP, takes forever to verify otherwise


%% collins_aerospace

disp("Running collins_aerospace..")

cab_path = vnncomp_path + "collins_aerospace_benchmark/";

cab_instances = ["onnx/yolov5nano_LRelu_640.onnx", "vnnlib/img_10012_perturbed_bbox_0_delta_0.001.vnnlib"];

% Run verification for collins_rul
for i=1:height(cab_instances)
    onnx = cab_path + cab_instances(i,1);
    vnnlib = cab_path + cab_instances(i,2);
    try
        run_vnncomp_instance("collins_aerospace_benchmark",onnx,vnnlib,"collins_aero_results_" + string(i)+".txt");
    catch ME
        warning("Failed")
        disp(onnx+"___"+vnnlib)
        warning(ME.message)
    end
end

% This instance finished (counterexample)


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
    try
        run_vnncomp_instance("collins_rul",onnx,vnnlib,"collins_rul_results_" + string(i)+".txt");
    catch ME
        warning("Failed")
        disp(onnx+"___"+vnnlib)
        warning(ME.message)
    end
end

% All good, see submission 234


%% Cora

disp("Running Cora benchmark..")

cora_path = vnncomp_path + "cora_2024/";

% We can probably do exact analysis from the beginning on all instances here
cora_instances = ["onnx/mnist-point.onnx", "vnnlib/mnist-img0.vnnlib";...
    "onnx/mnist-trades.onnx", "vnnlib/mnist-img1.vnnlib";...
    "onnx/mnist-set.onnx", "vnnlib/mnist-img20.vnnlib";...
    "onnx/svhn-point.onnx", "vnnlib/svhn-img84.vnnlib";...
    "onnx/svhn-set.onnx", "vnnlib/svhn-img273.vnnlib";...
    "onnx/svhn-trades.onnx", "vnnlib/svhn-img410.vnnlib";...
    "onnx/cifar10-set.onnx", "vnnlib/cifar10-img347.vnnlib";...
    "onnx/cifar10-point.onnx", "vnnlib/cifar10-img353.vnnlib";...
    "onnx/cifar10-trades.onnx", "vnnlib/cifar10-img423.vnnlib"...
    ];

% Run verification for CORA
for i=1:length(cora_instances)
    onnx = cora_path + cora_instances(i,1);
    vnnlib = cora_path + cora_instances(i,2);
    try
        run_vnncomp_instance("cora",onnx,vnnlib,"cora_results_" + string(i)+".txt");
    catch ME
        warning("Failed")
        disp(onnx+"___"+vnnlib)
        warning(ME.message)
    end
end


% All good


%% dist_shift

disp("Running dist_shift...")

dist_path = vnncomp_path + "dist_shift_2023/";

dist_instances = ["onnx/mnist_concat.onnx" ,"vnnlib/index2101_delta0.13.vnnlib";...
    "onnx/mnist_concat.onnx" ,"vnnlib/index9181_delta0.13.vnnlib";...
    "onnx/mnist_concat.onnx" ,"vnnlib/index4260_delta0.13.vnnlib";...
    ];

% Run verification for dist_shift
for i=1:length(dist_instances)
    onnx = dist_path + dist_instances(i,1);
    vnnlib = dist_path + dist_instances(i,2);
    try
        run_vnncomp_instance("dist_shift",onnx,vnnlib,"dist_results_" + string(i)+".txt");
    catch ME
        warning("Failed")
        disp(onnx+"___"+vnnlib)
        warning(ME.message)
    end
end

% No errors


%% linearizeNN

disp("Running linearizeNN...")

lin_path = vnncomp_path + "linearizenn_2024/";

lin_instances = ["onnx/AllInOne_10_10.onnx" ,"vnnlib/prop_10_10.vnnlib";...
    "onnx/AllInOne_30_30.onnx" ,"vnnlib/prop_30_30_0.vnnlib";...
    "onnx/AllInOne_50_50.onnx" ,"vnnlib/prop_50_50_0.vnnlib";...
    "onnx/AllInOne_50_120.onnx", "vnnlib/prop_50_120.vnnlib";...
    "onnx/AllInOne_80_30.onnx", "vnnlib/prop_80_30_3.vnnlib";...
    "onnx/AllInOne_120_120.onnx", "vnnlib/prop_120_120_0.vnnlib";...
    ]; 

% Run verification for dist_shift
for i=1:length(lin_instances)
    onnx = lin_path + lin_instances(i,1);
    vnnlib = lin_path + lin_instances(i,2);
    try
        run_vnncomp_instance("linearizenn",onnx,vnnlib,"linear_results_" + string(i)+".txt");
    catch ME
        warning("Failed")
        disp(onnx+"___"+vnnlib)
        warning(ME.message)
    end
end

% Looks like they all run, just need to do the CP method



%% lsnc_relu

disp("Running lsnc_relu...")

lsnc_path = vnncomp_path + "lsnc_relu/";

lsnc_instances = ["onnx/relu_quadrotor2d_state.onnx" ,"vnnlib/quadrotor2d_state_0.vnnlib";...
    ]; 

% Run verification for dist_shift
for i=1:height(lsnc_instances)
    onnx = lsnc_path + lsnc_instances(i,1);
    vnnlib = lsnc_path + lsnc_instances(i,2);
    try
        run_vnncomp_instance("lsnc_relu",onnx,vnnlib,"lsnc_results_" + string(i)+".txt");
    catch ME
        warning("Failed")
        disp(onnx+"___"+vnnlib)
        warning(ME.message)
    end
end

% Looks like they all run, just need to do the CP method
% Test to make sure we are dealing with properties correctly


%% malbeware

disp("Running malbeware...")

malbeware_path = vnncomp_path + "malbeware/";

malbeware_instances = ["onnx/malware_malimg_family_scaled_linear-25.onnx" ,"vnnlib/malbeware_family-Autorun.K_label-5_eps-2_idx-29.vnnlib";...
    "onnx/malware_malimg_family_scaled_linear-25.onnx", "vnnlib/malbeware_family-Swizzor.gen!E_label-20_eps-2_idx-103.vnnlib";...
    "onnx/malware_malimg_family_scaled_4-25.onnx", "vnnlib/malbeware_family-Autorun.K_label-5_eps-1_idx-25.vnnlib";...
]; 

% Run verification for dist_shift
for i=1:length(malbeware_instances)
    onnx = malbeware_path + malbeware_instances(i,1);
    vnnlib = malbeware_path + malbeware_instances(i,2);
    try
        run_vnncomp_instance("malbeware",onnx,vnnlib,"malbeware_results_" + string(i)+".txt");
    catch ME
        warning("Failed")
        disp(onnx+"___"+vnnlib)
        warning(ME.message)
    end
end


% too many sats, let's check the reshapes (needReshape = 2 should be fine)


%% metaroom

disp("Running metaroom..")

metaroom_path = vnncomp_path + "metaroom_2023/";

% We can probably do exact analysis from the beginning on all instances here
metaroom_instances = ["./onnx/6cnn_ry_92_15_no_custom_OP.onnx", "./vnnlib/spec_idx_147_eps_0.00000436.vnnlib";...
    "./onnx/4cnn_ry_47_7_no_custom_OP.onnx", "./vnnlib/spec_idx_23_eps_0.00000436.vnnlib";...
    "./onnx/4cnn_tz_94_15_no_custom_OP.onnx", "./vnnlib/spec_idx_90_eps_0.00001000.vnnlib";...
    ];

% Run verification for metaroom
for i=1:length(metaroom_instances)
    onnx = metaroom_path + metaroom_instances(i,1);
    vnnlib = metaroom_path + metaroom_instances(i,2);
    try
        run_vnncomp_instance("metaroom",onnx,vnnlib,"metaroom_results_" + string(i)+".txt");
    catch ME
        warning("Failed")
        disp(onnx+"___"+vnnlib)
        warning(ME.message)
    end
end



%% ml4acopf

% Skip for now, cannot initialize network

% disp("Running ml4acopf...")

% ml4acopf_path = vnncomp_path + "ml4acopf_2024/";
% 
% ml4acopf_instances = [...
%     % "onnx/14_ieee_ml4acopf-linear-residual.onnx" ,"vnnlib/14_ieee_prop9.vnnlib";...
%     "onnx/118_ieee_ml4acopf.onnx", "vnnlib/118_ieee_prop2.vnnlib";...
%     "onnx/300_ieee_ml4acopf-linear-nonresidual.onnx", "vnnlib/300_ieee_prop3.vnnlib";...
%     "onnx/14_ieee_ml4acopf-linear-nonresidual.onnx", "vnnlib/14_ieee_prop9.vnnlib";...
%     "onnx/14_ieee_ml4acopf.onnx","vnnlib/14_ieee_prop3.vnnlib";...
%     "onnx/118_ieee_ml4acopf-linear-nonresidual.onnx", "vnnlib/118_ieee_prop6.vnnlib";...
%     "onnx/300_ieee_ml4acopf.onnx", "vnnlib/300_ieee_prop3.vnnlib";...
% ]; 
% 
% % Run verification for dist_shift
% for i=1:length(ml4acopf_instances)
%     onnx = ml4acopf_path + ml4acopf_instances(i,1);
%     vnnlib = ml4acopf_path + ml4acopf_instances(i,2);
%     try
%         run_vnncomp_instance("ml4acopf",onnx,vnnlib,"ml4acopf_results_" + string(i)+".txt");
%     catch ME
%         warning("Failed")
%         disp(onnx+"___"+vnnlib)
%         warning(ME.message)
%     end
% end


%% nn4sys

disp("Running nn4sys...")

nn4sys_path = vnncomp_path + "nn4sys/";

nn4sys_instances = [... % all other networks are not supported...
    % "onnx/lindex.onnx","vnnlib/lindex_1.vnnlib";...
    % "onnx/lindex_deep.onnx", "vnnlib/lindex_200.vnnlib";...
    "onnx/pensieve_big_parallel.onnx", "vnnlib/pensieve_parallel_0.vnnlib";...
    "onnx/pensieve_small_simple.onnx", "vnnlib/pensieve_simple_1.vnnlib";...
    "onnx/pensieve_small_parallel.onnx", "vnnlib/pensieve_parallel_4.vnnlib";...
    "onnx/mscn_128d.onnx", "vnnlib/cardinality_0_1_128.vnnlib";...
    "onnx/mscn_128d_dual.onnx", "vnnlib/cardinality_1_1_128_dual.vnnlib";...
    "onnx/mscn_2048d.onnx", "vnnlib/cardinality_0_1_2048.vnnlib";...
    "onnx/mscn_2048d_dual.onnx", "vnnlib/cardinality_1_1_2048_dual.vnnlib";...
    ];

% Run verification for nn4sys 
for i=1:length(nn4sys_instances)
    onnx = nn4sys_path + nn4sys_instances(i,1);
    vnnlib = nn4sys_path + nn4sys_instances(i,2);
    try
        run_vnncomp_instance("nn4sys",onnx,vnnlib,"nn4sys_results_" + string(i)+".txt");
    catch ME
        warning("Failed")
        disp(onnx+"___"+vnnlib)
        warning(ME.message)
    end
end

% Everything looks okey, we need to test cp now


%% relusplitter

disp("Running relusplitter...")

relusplitter_path = vnncomp_path + "relusplitter/";

relusplitter_instances = [... 
    "onnx/mnist_fc_vnncomp2022_mnist-net_256x4.onnx","vnnlib/mnist_fc_vnncomp2022_prop_7_0.05.vnnlib";...
    "onnx/mnist_fc_vnncomp2022_mnist-net_256x4_RSPLITTER_mnist_fc_vnncomp2022_prop_7_0.05.onnx", "vnnlib/mnist_fc_vnncomp2022_prop_7_0.05.vnnlib";...
    "onnx/mnist_fc_vnncomp2022_mnist-net_256x6.onnx", "vnnlib/mnist_fc_vnncomp2022_prop_4_0.03.vnnlib";...
    "onnx/mnist_fc_vnncomp2022_mnist-net_256x6_RSPLITTER_mnist_fc_vnncomp2022_prop_5_0.05.onnx", "vnnlib/mnist_fc_vnncomp2022_prop_6_0.03.vnnlib";...
    "onnx/oval21-benchmark_cifar_base_kw.onnx", "vnnlib/oval21-benchmark_cifar_base_kw-img2069-eps0.02718954248366013.vnnlib";...
    "onnx/oval21-benchmark_cifar_base_kw_RSPLITTER_oval21-benchmark_cifar_base_kw-img5541-eps0.010718954248366015.onnx", "vnnlib/oval21-benchmark_cifar_base_kw-img5541-eps0.010718954248366015.vnnlib";...
    "onnx/oval21-benchmark_cifar_deep_kw.onnx", "vnnlib/oval21-benchmark_cifar_deep_kw-img221-eps0.025098039215686277.vnnlib";...
    "onnx/oval21-benchmark_cifar_deep_kw_RSPLITTER_oval21-benchmark_cifar_deep_kw-img6508-eps0.044183006535947714.onnx", "vnnlib/oval21-benchmark_cifar_deep_kw-img6508-eps0.044183006535947714.vnnlib";...
    "onnx/cifar_biasfield_vnncomp2022_cifar_bias_field_5.onnx", "vnnlib/cifar_biasfield_vnncomp2022_prop_5.vnnlib";...
    "onnx/cifar_biasfield_vnncomp2022_cifar_bias_field_27_RSPLITTER_cifar_biasfield_vnncomp2022_prop_27.onnx", "vnnlib/cifar_biasfield_vnncomp2022_prop_27.vnnlib";...
    ];

% Run verification for relusplitter
for i=1:length(relusplitter_instances)
    onnx = relusplitter_path + relusplitter_instances(i,1);
    vnnlib = relusplitter_path + relusplitter_instances(i,2);
    try
        run_vnncomp_instance("relusplitter",onnx,vnnlib,"relusplitter_results_" + string(i)+".txt");
    catch ME
        warning("Failed")
        disp(onnx+"___"+vnnlib)
        warning(ME.message)
    end
end

% oval seem fine
% bias field super slow, let's use cp for that


%% safeNLP

disp("Running safeNLP..")

safeNLP_path = vnncomp_path + "safenlp_2024/";

% We can probably do exact analysis from the beginning on all instances here
safeNLP_instances = ["onnx/ruarobot/perturbations_0.onnx", "vnnlib/ruarobot/hyperrectangle_3607.vnnlib";...
    "onnx/medical/perturbations_0.onnx", "vnnlib/medical/hyperrectangle_369.vnnlib";...
    ];

% Run verification for safeNLP
for i=1:length(safeNLP_instances)
    onnx = safeNLP_path + safeNLP_instances(i,1);
    vnnlib = safeNLP_path + safeNLP_instances(i,2);
    try
        run_vnncomp_instance("safenlp",onnx,vnnlib,"safeNLP_results_" + string(i)+".txt");
    catch ME
        warning("Failed")
        disp(onnx+"___"+vnnlib)
        warning(ME.message)
    end
end

% very short timeouts, barely any time to compute anything
% current approach is okay
% let's try cp to see if any better


%% sat relu

disp("Running sat_relu..")

satrelu_path = vnncomp_path + "sat_relu/";

% We can probably do exact analysis from the beginning on all instances here
satrelu_instances = ["onnx/sat_v18_c75.onnx", "vnnlib/sat_v18_c75.vnnlib";...
    "onnx/unsat_v18_c75.onnx", "vnnlib/unsat_v18_c75.vnnlib";...
    ];

% Run verification for satrelu
for i=1:length(satrelu_instances)
    onnx = satrelu_path + satrelu_instances(i,1);
    vnnlib = satrelu_path + satrelu_instances(i,2);
    try
        run_vnncomp_instance("sat_relu",onnx,vnnlib,"satrelu_results_" + string(i)+".txt");
    catch ME
        warning("Failed")
        disp(onnx+"___"+vnnlib)
        warning(ME.message)
    end
end

% very fast computation with relax, but not many verified
% Let's test on submission site


%% soundnessbench

disp("Running soundness benchmark..")

sound_path = vnncomp_path + "soundnessbench/";

% We can probably do exact analysis from the beginning on all instances here
sound_instances = ["onnx/model.onnx", "vnnlib/model_0.vnnlib";...
    "onnx/model.onnx", "vnnlib/model_29.vnnlib";...
    ];

% Run verification for satrelu
for i=1:length(sound_instances)
    onnx = sound_path + sound_instances(i,1);
    vnnlib = sound_path + sound_instances(i,2);
    try
        run_vnncomp_instance("soundnessbench",onnx,vnnlib,"sound_results_" + string(i)+".txt");
    catch ME
        warning("Failed")
        disp(onnx+"___"+vnnlib)
        warning(ME.message)
    end
end

% This one runs out of memory even with the relax star method for all
% instances (weird...)


%% tinyimagenet

disp("Running tinyimagenet..")

tinyimagenet_path = vnncomp_path + "tinyimagenet_2024/";

% no vnnlib properties provided in this year's repo
% (Some may be available in 2022 repo)
tinyimagenet_instances = ["onnx/TinyImageNet_resnet_medium.onnx", "vnnlib/TinyImageNet_resnet_medium_prop_idx_762_sidx_4897_eps_0.0039.vnnlib"; ...
    ];

% Run verification for tinyimagenet
for i=1:height(tinyimagenet_instances)
    onnx = tinyimagenet_path + tinyimagenet_instances(i,1);
    vnnlib = tinyimagenet_path + tinyimagenet_instances(i,2);
    try
        run_vnncomp_instance("tinyimagenet",onnx,vnnlib,"tinyimagenet_results_" + string(i)+".txt");
    catch ME
        warning("Failed")
        disp(onnx+"___"+vnnlib)
        warning(ME.message)
    end
end

% out of memory


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
    try
        run_vnncomp_instance("tllverifybench",onnx,vnnlib,"tllverify_results_" + string(i)+".txt");
    catch ME
        warning("Failed")
        disp(onnx+"___"+vnnlib)
        warning(ME.message)
    end
end

% all unknown or sat

%% traffic signs recognition

disp("Running traffic signs recognition..")

traffic_path = vnncomp_path + "traffic_signs_recognition_2023/";

traffic_instances = ["onnx/3_30_30_QConv_16_3_QConv_32_2_Dense_43_ep_30.onnx", "vnnlib/model_30_idx_7377_eps_1.00000.vnnlib";...
    "onnx/3_48_48_QConv_32_5_MP_2_BN_QConv_64_5_MP_2_BN_QConv_64_3_BN_Dense_256_BN_Dense_43_ep_30.onnx", "vnnlib/model_48_idx_7377_eps_1.00000.vnnlib";...
    "onnx/3_64_64_QConv_32_5_MP_2_BN_QConv_64_5_MP_2_BN_QConv_64_3_MP_2_BN_Dense_1024_BN_Dense_43_ep_30.onnx","vnnlib/model_64_idx_7377_eps_1.00000.vnnlib";...
    ]; 

% Run verification for traffic
for i=1:length(traffic_instances)
    onnx = traffic_path + traffic_instances(i,1);
    vnnlib = traffic_path + traffic_instances(i,2);
    try
        run_vnncomp_instance("traffic_signs_recognition",onnx,vnnlib,"traffic_results_" + string(i)+".txt");
    catch ME
        warning("Failed")
        disp(onnx+"___"+vnnlib)
        warning(ME.message)
    end
end

% Fine, let's test how fast we can do CP on submission site
% We got this warning tho: 
% Warning: Data format "BCSS" you specified for input 1 does not match the format "BSSC" derived by the software. The data format derived
% by the software will be used. 

%% vggnet

disp("Running vggnet16...")

vgg_path = vnncomp_path + "vggnet16_2022/";

vgg_instances = ["onnx/vgg16-7.onnx", "vnnlib/spec0_briard.vnnlib"]; 

% Run verification for vgg
for i=1:height(vgg_instances)
    onnx = vgg_path + vgg_instances(i,1);
    vnnlib = vgg_path + vgg_instances(i,2);
    try
        run_vnncomp_instance("vggnet16",onnx,vnnlib,"vgg_results_" + string(i)+".txt");
    catch ME
        warning("Failed")
        disp(onnx+"___"+vnnlib)
        warning(ME.message)
    end
end


%% vit

disp("Running vit...")

vit_path = vnncomp_path + "vit_2023/";

vit_instances = ["onnx/pgd_2_3_16.onnx", "vnnlib/pgd_2_3_16_8835.vnnlib";...
    "onnx/ibp_3_3_8.onnx", "vnnlib/ibp_3_3_8_4868.vnnlib";...
    ]; 

% Run verification for vit
for i=1:length(vit_instances)
    onnx = vit_path + vit_instances(i,1);
    vnnlib = vit_path + vit_instances(i,2);
    try
        run_vnncomp_instance("vit_2023",onnx,vnnlib,"vit_results_" + string(i)+".txt");
    catch ME
        warning("Failed")
        disp(onnx+"___"+vnnlib)
        warning(ME.message)
    end
end

%% yolo

disp("Running yolo...")

yolo_path = vnncomp_path + "yolo_2023/";

yolo_instances = ["onnx/TinyYOLO.onnx", "vnnlib/TinyYOLO_prop_000034_eps_1_255.vnnlib"]; 

% Run verification for yolo
for i=1:height(yolo_instances)
    onnx = yolo_path + yolo_instances(i,1);
    vnnlib = yolo_path + yolo_instances(i,2);
    try
        run_vnncomp_instance("yolo_2023",onnx,vnnlib,"yolo_results_" + string(i)+".txt");
    catch ME
        warning("Failed")
        disp(onnx+"___"+vnnlib)
        warning(ME.message)
    end
end