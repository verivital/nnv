function Run_iteration(i)


%% tiny image

disp("Running tinyimagenet..")



dir =  pwd;
addpath(genpath(dir));

parentDir = fileparts(dir);
src_dir = parentDir + "\src";
addpath(genpath(src_dir))

tinyimagenet_instances = instance_generator("");

onnx = tinyimagenet_instances(i,1);
vnnlib = tinyimagenet_instances(i,2);
run_vnncomp2025("tinyimagenet",onnx,vnnlib,"tinyimagenet_results_" + string(i)+".txt");


rmpath(genpath(dir))
end

%% It takes a long time to perform reachability