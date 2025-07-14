function Run_iteration(i)


%% Sound

disp("Running Sound..")

dir =  pwd;
addpath(genpath(dir));


parentDir = fileparts(dir);
src_dir = parentDir + "\src";
addpath(genpath(src_dir))

% Sound_instances = instance_generator("");

Sound_instances = {
    "onnx\model.onnx", "vnnlib\model_0.vnnlib";
    "onnx\model.onnx", "vnnlib\model_1.vnnlib";
    "onnx\model.onnx", "vnnlib\model_2.vnnlib";
    "onnx\model.onnx", "vnnlib\model_3.vnnlib";
    "onnx\model.onnx", "vnnlib\model_4.vnnlib";
    "onnx\model.onnx", "vnnlib\model_5.vnnlib";
    "onnx\model.onnx", "vnnlib\model_6.vnnlib";
    "onnx\model.onnx", "vnnlib\model_7.vnnlib";
    "onnx\model.onnx", "vnnlib\model_8.vnnlib";
    "onnx\model.onnx", "vnnlib\model_9.vnnlib";
    "onnx\model.onnx", "vnnlib\model_10.vnnlib";
    "onnx\model.onnx", "vnnlib\model_11.vnnlib";
    "onnx\model.onnx", "vnnlib\model_12.vnnlib";
    "onnx\model.onnx", "vnnlib\model_13.vnnlib";
    "onnx\model.onnx", "vnnlib\model_14.vnnlib";
    "onnx\model.onnx", "vnnlib\model_15.vnnlib";
    "onnx\model.onnx", "vnnlib\model_16.vnnlib";
    "onnx\model.onnx", "vnnlib\model_17.vnnlib";
    "onnx\model.onnx", "vnnlib\model_18.vnnlib";
    "onnx\model.onnx", "vnnlib\model_19.vnnlib";
    "onnx\model.onnx", "vnnlib\model_20.vnnlib";
    "onnx\model.onnx", "vnnlib\model_21.vnnlib";
    "onnx\model.onnx", "vnnlib\model_22.vnnlib";
    "onnx\model.onnx", "vnnlib\model_23.vnnlib";
    "onnx\model.onnx", "vnnlib\model_24.vnnlib";
    "onnx\model.onnx", "vnnlib\model_25.vnnlib";
    "onnx\model.onnx", "vnnlib\model_26.vnnlib";
    "onnx\model.onnx", "vnnlib\model_27.vnnlib";
    "onnx\model.onnx", "vnnlib\model_28.vnnlib";
    "onnx\model.onnx", "vnnlib\model_29.vnnlib";
    "onnx\model.onnx", "vnnlib\model_30.vnnlib";
    "onnx\model.onnx", "vnnlib\model_31.vnnlib";
    "onnx\model.onnx", "vnnlib\model_32.vnnlib";
    "onnx\model.onnx", "vnnlib\model_33.vnnlib";
    "onnx\model.onnx", "vnnlib\model_34.vnnlib";
    "onnx\model.onnx", "vnnlib\model_35.vnnlib";
    "onnx\model.onnx", "vnnlib\model_36.vnnlib";
    "onnx\model.onnx", "vnnlib\model_37.vnnlib";
    "onnx\model.onnx", "vnnlib\model_38.vnnlib";
    "onnx\model.onnx", "vnnlib\model_39.vnnlib";
    "onnx\model.onnx", "vnnlib\model_40.vnnlib";
    "onnx\model.onnx", "vnnlib\model_41.vnnlib";
    "onnx\model.onnx", "vnnlib\model_42.vnnlib";
    "onnx\model.onnx", "vnnlib\model_43.vnnlib";
    "onnx\model.onnx", "vnnlib\model_44.vnnlib";
    "onnx\model.onnx", "vnnlib\model_45.vnnlib";
    "onnx\model.onnx", "vnnlib\model_46.vnnlib";
    "onnx\model.onnx", "vnnlib\model_47.vnnlib";
    "onnx\model.onnx", "vnnlib\model_48.vnnlib";
    "onnx\model.onnx", "vnnlib\model_49.vnnlib";
};

onnx = Sound_instances(i,1);
onnx = onnx{1};
vnnlib = Sound_instances(i,2);
vnnlib = vnnlib{1};
run_vnncomp2025("Sound",onnx,vnnlib,"Sound_results_" + string(i)+".txt");

rmpath(genpath(dir))
end


%% It takes for a while to finish reachability