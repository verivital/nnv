%% Run some verification instances from the competition

% path to benchmarks
% vnncomp_path = "/home/manzand/Documents/MATLAB/vnncomp2023_benchmarks/benchmarks/";
vnncomp_path = "/home/dieman95/Documents/MATLAB/vnncomp2023_benchmarks/benchmarks/";

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

% Run verification for acas
for i=1:length(vit_instances)
    onnx = vit_path + vit_instances(i,1);
    vnnlib = vit_path + vit_instances(i,2);
    run_vnncomp_instance("vit",onnx,vnnlib,'results.txt');
end


%% nn4sys


%% dist_shift


%% traffic_sign
% we should be able to support it, but matlab does not, so we would have to create the models manually...
% can load most info with Keras importer

%% collins_rul


%% cgan


%% vggnet


%% ml4acopf


%% tllverify

