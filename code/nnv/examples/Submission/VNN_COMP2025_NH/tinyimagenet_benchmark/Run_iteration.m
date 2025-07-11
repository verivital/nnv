function Run_iteration(i)


%% Sound

disp("Running tinyimagenet..")


dir =  pwd;
parentDir = fileparts(dir);
src_dir = parentDir + "\src";
addpath(genpath(src_dir))

parentDir = fileparts(parentDir);
parentDir = fileparts(parentDir);
nnv_dir = fileparts(parentDir);
addpath(genpath(nnv_dir))

tinyimagenet_instances = instance_generator("");

onnx = tinyimagenet_instances(i,1);
vnnlib = tinyimagenet_instances(i,2);
run_vnncomp2025("tinyimagenet",onnx,vnnlib,"tinyimagenet_results_" + string(i)+".txt");


end

%% It takes a long time to perform reachability