function Run_iteration(i)


%% LSNC

disp("Running LCNS_ReLU..")

dir =  pwd;
addpath(genpath(dir));


parentDir = fileparts(dir);
src_dir = parentDir + "\src";
addpath(genpath(src_dir))

LCNS_instances = instance_generator("");

onnx = LCNS_instances(i,1);
vnnlib = LCNS_instances(i,2);
run_vnncomp2025("LCNS_ReLU",onnx,vnnlib,"LCNS_results_" + string(i)+".txt");

rmpath(genpath(dir))

end
