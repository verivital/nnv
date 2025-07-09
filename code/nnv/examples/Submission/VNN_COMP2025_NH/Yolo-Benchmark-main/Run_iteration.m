function Run_iteration(i)


%% yolo

disp("Running yolo..")


dir =  pwd;
parentDir = fileparts(dir);
src_dir = parentDir + "\src";
addpath(genpath(src_dir))

parentDir = fileparts(parentDir);
parentDir = fileparts(parentDir);
nnv_dir = fileparts(parentDir);
addpath(genpath(nnv_dir))

yolo_instances = instance_generator("");

onnx = yolo_instances(i,1);
vnnlib = yolo_instances(i,2);
run_vnncomp2025("yolo",onnx,vnnlib,"yolo_results_" + string(i)+".txt");


end


%% It takes a long time to perform reachability