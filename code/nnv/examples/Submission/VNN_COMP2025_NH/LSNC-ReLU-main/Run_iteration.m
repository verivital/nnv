function Run_iteration(i)


%% cgan

disp("Running LCNS_ReLU..")

dir =  pwd;
parentDir = fileparts(dir);
src_dir = parentDir + "\src";
addpath(genpath(src_dir))


parentDir = fileparts(parentDir);
parentDir = fileparts(parentDir);
nnv_dir = fileparts(parentDir);
addpath(genpath(nnv_dir))


LCNS_instances = instance_generator("");

onnx = LCNS_instances(i,1);
vnnlib = LCNS_instances(i,2);
run_vnncomp2025("LCNS_ReLU",onnx,vnnlib,"LCNS_results_" + string(i)+".txt");


end



%% Layer 8 is a relu_quadrotor2d_state.MatMul_To_AddLayer1017 which have not supported yet in nnv, please consider removing this layer for the analysis --->  MatMulToAddLayer
%% Layer 11 is a relu_quadrotor2d_state.MatMul_To_ReduceSumLayer1005 ---> MatMulToReduceSumLayer
%% Layer 12 is a relu_quadrotor2d_state.MatMul_To_SubLayer1025 ---> 


%% The MatMul Layer was not clear what is the underlying operation. The onnx format was entirely not helpful in Netron sice the dlnetwork version was so different , on the other
%% the dlnetwork version did not give us the underlying formula as I was not able to open the details in the analyzeNetwork()%%

