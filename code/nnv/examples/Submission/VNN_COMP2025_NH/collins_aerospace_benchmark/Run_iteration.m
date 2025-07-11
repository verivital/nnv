function Run_iteration(i)


%% cgan

disp("Running collins..")

dir =  pwd;
parentDir = fileparts(dir);
src_dir = parentDir + "\src";
addpath(genpath(src_dir))


parentDir = fileparts(parentDir);
parentDir = fileparts(parentDir);
nnv_dir = fileparts(parentDir);
addpath(genpath(nnv_dir))

% cgan_instances = instance_generator("");
cgan_instances = {
    "onnx/yolov5nano_LRelu_640.onnx",	"vnnlib/img_10012_perturbed_bbox_0_delta_0.001.vnnlib";
    "onnx/yolov5nano_LRelu_640.onnx", 	"vnnlib/img_398_perturbed_bbox_4_delta_0.005.vnnlib";
    "onnx/yolov5nano_LRelu_640.onnx",	"vnnlib/img_15588_perturbed_bbox_0_delta_0.01.vnnlib";
    "onnx/yolov5nano_LRelu_640.onnx", 	"vnnlib/img_14147_perturbed_bbox_2_delta_0.02.vnnlib";
    "onnx/yolov5nano_LRelu_640.onnx",	"vnnlib/img_14147_perturbed_bbox_1_delta_0.05.vnnlib";
    "onnx/yolov5nano_LRelu_640.onnx",	"vnnlib/img_6876_perturbed_bbox_3_delta_0.1.vnnlib";
};

onnx = cgan_instances(i,1);
onnx = onnx{1};
vnnlib = cgan_instances(i,2);
vnnlib = vnnlib{1};
run_vnncomp2025("collins_aerospace_benchmark",onnx,vnnlib,"collins_results_" + string(i)+".txt");


end

%% Layer 86 is a nnet.cnn.layer.Resize2DLayer which have not supported yet in nnv, please consider removing this layer for the analysis --> Resize2DLayer
%% Layer 145 is a yolov5nano_LRelu_640.Reshape_To_ConcatLayer1159 -->  ReshapeToConcatenationLayer