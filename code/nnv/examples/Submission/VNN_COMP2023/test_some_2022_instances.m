%% Run some verification instances from the competition

% path to benchmarks
% vnncomp_path = "/home/manzand/Documents/MATLAB/vnncomp2022_benchmarks/benchmarks/";
vnncomp_path = "/home/dieman95/Documents/MATLAB/vnncomp2022_benchmarks/benchmarks/";

% Go through some of the instances for every benchmark


%% MNIST fc

disp("Running mnist_fc")

bench_path = vnncomp_path + "mnist_fc/";

bench_instances = [
    "onnx/mnist-net_256x2.onnx" , "vnnlib/prop_14_0.05.vnnlib";...
    "onnx/mnist-net_256x4.onnx" , "vnnlib/prop_0_0.03.vnnlib";...
%     "onnx/mnist-net_256x6.onnx" , "vnnlib/prop_14_0.05.vnnlib";... % takes too long this one
]; 

% Run verification for acas 
for i=1:length(bench_instances)
    onnx = bench_path + bench_instances(i,1);
    vnnlib = bench_path + bench_instances(i,2);
    try
        run_vnncomp2023_instance("mnist_fc",onnx,vnnlib,"mnist_results_" + string(i)+".txt");
    end
end


%% cifar100_tinyimagenet_resnet

disp("Running cifar100 tinyimagenet...")

bench_path = vnncomp_path + "cifar100_tinyimagenet_resnet/";

bench_instances = [
    "onnx/CIFAR100_resnet_small.onnx" , "vnnlib/CIFAR100_resnet_small_prop_idx_1899_sidx_9565_eps_0.0039.vnnlib";...
    "onnx/CIFAR100_resnet_medium.onnx" , "vnnlib/CIFAR100_resnet_medium_prop_idx_7764_sidx_254_eps_0.0039.vnnlib";...
    "onnx/CIFAR100_resnet_super.onnx" , "vnnlib/CIFAR100_resnet_super_prop_idx_9014_sidx_9195_eps_0.0039.vnnlib";...
    "onnx/TinyImageNet_resnet_medium.onnx" , "vnnlib/TinyImageNet_resnet_medium_prop_idx_3688_sidx_720_eps_0.0039.vnnlib"]; 

% Run verification for acas 
for i=1:length(bench_instances)
    onnx = bench_path + bench_instances(i,1);
    vnnlib = bench_path + bench_instances(i,2);
    try
        run_vnncomp2023_instance("cifar100_tinyimagenet_resnet",onnx,vnnlib,"cifar100_results_" + string(i)+".txt");
    end
end

%% cifar2020

disp("Running cifar2020...")

bench_path = vnncomp_path + "cifar2020/";

bench_instances = [
    "onnx/cifar10_2_255_simplified.onnx", "vnnlib/cifar10_spec_idx_99_eps_0.00784_n1.vnnlib";...
    "onnx/cifar10_8_255_simplified.onnx" , "vnnlib/cifar10_spec_idx_1_eps_0.03137_n1.vnnlib";...
    "onnx/convBigRELU__PGD.onnx" , "vnnlib/cifar10_spec_idx_1_eps_0.00784.vnnlib"]; 

% Run verification for acas 
for i=1:length(bench_instances)
    onnx = bench_path + bench_instances(i,1);
    vnnlib = bench_path + bench_instances(i,2);
    try
        run_vnncomp2023_instance("cifar2020",onnx,vnnlib,"cifar2020_results_" + string(i)+".txt");
    end
end

%% cifar biasfield

disp("Running cifar biasfield...")

bench_path = vnncomp_path + "cifar_biasfield/";

bench_instances = [
    "onnx/cifar_bias_field_0.onnx" , "vnnlib/prop_0.vnnlib";...
    "onnx/cifar_bias_field_1.onnx" , "vnnlib/prop_1.vnnlib";...
    "onnx/cifar_bias_field_2.onnx" , "vnnlib/prop_2.vnnlib"]; 

% Run verification for acas 
for i=1:length(bench_instances)
    onnx = bench_path + bench_instances(i,1);
    vnnlib = bench_path + bench_instances(i,2);
    try
        run_vnncomp2023_instance("cifar_biasfield",onnx,vnnlib,"cifarBiasfield_results_" + string(i)+".txt");
    end
end

%% oval21

disp("Running oval21..")

bench_path = vnncomp_path + "oval21/";

bench_instances = [
    "onnx/cifar_wide_kw.onnx" , "vnnlib/cifar_wide_kw-img3458-eps0.022875816993464054.vnnlib";...
    "onnx/cifar_deep_kw.onnx" , "vnnlib/cifar_deep_kw-img4405-eps0.036732026143790855.vnnlib";...
    "onnx/cifar_base_kw.onnx" , "vnnlib/cifar_base_kw-img8194-eps0.018300653594771243.vnnlib";...
    ]; 

% Run verification for acas 
for i=1:length(bench_instances)
    onnx = bench_path + bench_instances(i,1);
    vnnlib = bench_path + bench_instances(i,2);
    try
        run_vnncomp2023_instance("oval21",onnx,vnnlib,"oval21_results_" + string(i)+".txt");
    end
end

%% reach_prob_density

disp("Running reach_prob_density..")

bench_path = vnncomp_path + "reach_prob_density/";

bench_instances = [
    "onnx/vdp.onnx" , "vnnlib/vdp_11.vnnlib";...
    "onnx/robot.onnx" , "vnnlib/robot_0.vnnlib";...
    "onnx/gcas.onnx" , "vnnlib/gcas_0.vnnlib";...
    ]; 

% Run verification for acas 
for i=1:length(bench_instances)
    onnx = bench_path + bench_instances(i,1);
    vnnlib = bench_path + bench_instances(i,2);
    try
        run_vnncomp2023_instance("reach_prob_density",onnx,vnnlib,"reach_prob_results_" + string(i)+".txt");
    end
end

%% rl benchmarks

disp("Running RL benchmarks..")

bench_path = vnncomp_path + "rl_benchmarks/";

bench_instances = [...
    "onnx/cartpole.onnx" , "vnnlib/cartpole_case_unsafe_8.vnnlib";...
    "onnx/cartpole.onnx" , "vnnlib/cartpole_case_safe_9.vnnlib";...
    "onnx/lunarlander.onnx" , "vnnlib/lunarlander_case_safe_2.vnnlib";...
    "onnx/lunarlander.onnx" , "vnnlib/lunarlander_case_unsafe_3.vnnlib";...
    "onnx/dubinsrejoin.onnx" , "vnnlib/dubinsrejoin_case_safe_6.vnnlib";...
    "onnx/dubinsrejoin.onnx" , "vnnlib/dubinsrejoin_case_unsafe_7.vnnlib";...
    ]; 

% Run verification for acas 
for i=1:length(bench_instances)
    onnx = bench_path + bench_instances(i,1);
    vnnlib = bench_path + bench_instances(i,2);
    try
        run_vnncomp2023_instance("rl_benchmarks",onnx,vnnlib,"rl_results_" + string(i)+".txt");
    end
end


%% sri resnets

disp("Running SRI benchmarks..")

bench_path = vnncomp_path;

bench_instances = [...
    "sri_resnet_a/onnx/resnet_3b2_bn_mixup_adv_4.0_bs128_lr-1.onnx" , "sri_resnet_a/vnnlib/cifar10_spec_idx_1551_eps_0.00198.vnnlib";...
    "sri_resnet_b/onnx/resnet_3b2_bn_mixup_ssadv_4.0_bs128_lr-1_v2.onnx" , "sri_resnet_b/vnnlib/cifar10_spec_idx_3684_eps_0.00350.vnnlib";...
    ]; 

% Run verification for acas 
for i=1:length(bench_instances)
    onnx = bench_path + bench_instances(i,1);
    vnnlib = bench_path + bench_instances(i,2);
%     try
    run_vnncomp2023_instance("sri_resnet",onnx,vnnlib,"sri_results_" + string(i)+".txt");
%     end
end

disp("Finished running all test instances...");
