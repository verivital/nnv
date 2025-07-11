function Run_iteration(i)


%% ml4acopf

disp("Running ml4acopf..")

dir =  pwd;
parentDir = fileparts(dir);
src_dir = parentDir + "\src";
addpath(genpath(src_dir))


parentDir = fileparts(parentDir);
parentDir = fileparts(parentDir);
nnv_dir = fileparts(parentDir);
addpath(genpath(nnv_dir))

ml4acopf_instances = instance_generator("");

onnx = ml4acopf_instances(i,1);
vnnlib = ml4acopf_instances(i,2);
run_vnncomp2025("ml4acopf",onnx,vnnlib,"ML4aCOPF_results_" + string(i)+".txt");


end

%% Layer 10 is a x118_ieee_ml4acopf_linear_nonresidual.Slice_To_MatMulLayer1007 which have not supported yet in nnv, please consider removing this layer for the analysis 
