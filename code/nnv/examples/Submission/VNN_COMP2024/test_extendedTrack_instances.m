%% Run some verification instances from the competition

% path to benchmarks
vnncomp_path = "/home/manzand/Documents/MATLAB/vnncomp2023_benchmarks/benchmarks/"; % change to whatever local path vnncomp benchmarks are stored
% vnncomp_path = "/home/dieman95/Documents/MATLAB/vnncomp2023_benchmarks/benchmarks/";


%% ViT

disp("Running ViT...")

vit_path = vnncomp_path + "vit/";

vit_instances = [...
    "onnx/ibp_3_3_8.onnx" , "vnnlib/ibp_3_3_8_3850.vnnlib";...
    "onnx/pgd_2_3_16.onnx","vnnlib/pgd_2_3_16_1744.vnnlib";...
    ];

% Run verification for vit 
for i=1:length(vit_instances)
    onnx = vit_path + vit_instances(i,1);
    vnnlib = vit_path + vit_instances(i,2);
    run_vnncomp2024_instance("vit",onnx,vnnlib,"vit_results_" + string(i)+".txt");
end



%% traffic_sign
% we should be able to support it, but matlab does not, so we would have to create the models manually...
% can load most info with Keras importer


%% vggnet16

disp("Running vggnet16...")

vgg_path = vnncomp_path + "vggnet16/";

vgg_instances = ["onnx/vgg16-7.onnx","vnnlib/spec0_screw.vnnlib"];

% Run verification for vgg16 (error when creating random examples, array size too large)
% I can imagine that similar errors will be encountered when running verification
onnx = vgg_path + vgg_instances(1,1);
vnnlib = vgg_path + vgg_instances(1,2);
run_vnncomp2024_instance("vggnet16",onnx,vnnlib,"vgg_results.txt");


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

