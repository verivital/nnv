function Run_iteration(i)


%% collins

disp("Running collins..")
dir =  pwd;
addpath(genpath(dir));

parentDir = fileparts(dir);
src_dir = parentDir + "\src";
addpath(genpath(src_dir))

% cgan_instances = instance_generator("");

collins_instances = {
    "onnx/yolov5nano_LRelu_640.onnx",	"vnnlib/img_10012_perturbed_bbox_0_delta_0.001.vnnlib";
    "onnx/yolov5nano_LRelu_640.onnx", 	"vnnlib/img_398_perturbed_bbox_4_delta_0.005.vnnlib";
    "onnx/yolov5nano_LRelu_640.onnx",	"vnnlib/img_15588_perturbed_bbox_0_delta_0.01.vnnlib";
    "onnx/yolov5nano_LRelu_640.onnx", 	"vnnlib/img_14147_perturbed_bbox_2_delta_0.02.vnnlib";
    "onnx/yolov5nano_LRelu_640.onnx",	"vnnlib/img_14147_perturbed_bbox_1_delta_0.05.vnnlib";
    "onnx/yolov5nano_LRelu_640.onnx",	"vnnlib/img_6876_perturbed_bbox_3_delta_0.1.vnnlib";
};

onnx = collins_instances(i,1);
onnx = onnx{1};
vnnlib = collins_instances(i,2);
vnnlib = vnnlib{1};
run_vnncomp2025("collins_aerospace_benchmark",onnx,vnnlib,"collins_results_" + string(i)+".txt");

rmpath(genpath(dir))

end