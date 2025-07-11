function Run_iteration(i)


%% ViT

disp("Running ViT..")

dir =  pwd;
parentDir = fileparts(dir);
src_dir = parentDir + "\src";
addpath(genpath(src_dir))

parentDir = fileparts(parentDir);
parentDir = fileparts(parentDir);
nnv_dir = fileparts(parentDir);
addpath(genpath(nnv_dir))

ViT_instances = instance_generator("");

onnx = ViT_instances(i,1);
vnnlib = ViT_instances(i,2);
run_vnncomp2025("ViT",onnx,vnnlib,"yolo_results_" + string(i)+".txt");


end


%% Layer 3 is a ibp_3_3_8.Shape_To_ReduceMeanLayer1157 which have not supported yet in nnv, please consider removing this layer for the analysis 
