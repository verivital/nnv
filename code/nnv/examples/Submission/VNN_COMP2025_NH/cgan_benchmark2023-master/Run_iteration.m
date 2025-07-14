function Run_iteration(i)


%% cgan
addpath(genpath('C:\Users\navid\Documents\nnv\code'))
disp("Running cgan..")

dir =  pwd;
addpath(genpath(dir));
parentDir = fileparts(dir);
src_dir = parentDir + "\src";
addpath(genpath(src_dir))

cgan_instances = {
"onnx/cGAN_imgSz32_nCh_1.onnx",	"vnnlib/cGAN_imgSz32_nCh_1_prop_0_input_eps_0.010_output_eps_0.015.vnnlib";
"onnx/cGAN_imgSz32_nCh_1.onnx",	"vnnlib/cGAN_imgSz32_nCh_1_prop_1_input_eps_0.010_output_eps_0.015.vnnlib";
"onnx/cGAN_imgSz32_nCh_1.onnx",	"vnnlib/cGAN_imgSz32_nCh_1_prop_2_input_eps_0.020_output_eps_0.025.vnnlib";
"onnx/cGAN_imgSz32_nCh_1.onnx",	"vnnlib/cGAN_imgSz32_nCh_1_prop_3_input_eps_0.020_output_eps_0.025.vnnlib";
"onnx/cGAN_imgSz32_nCh_3.onnx",	"vnnlib/cGAN_imgSz32_nCh_3_prop_0_input_eps_0.015_output_eps_0.020.vnnlib";
"onnx/cGAN_imgSz32_nCh_3.onnx",	"vnnlib/cGAN_imgSz32_nCh_3_prop_1_input_eps_0.015_output_eps_0.020.vnnlib";
"onnx/cGAN_imgSz32_nCh_3.onnx",	"vnnlib/cGAN_imgSz32_nCh_3_prop_2_input_eps_0.010_output_eps_0.015.vnnlib";
"onnx/cGAN_imgSz32_nCh_3.onnx",	"vnnlib/cGAN_imgSz32_nCh_3_prop_3_input_eps_0.015_output_eps_0.020.vnnlib";
"onnx/cGAN_imgSz64_nCh_1.onnx",	"vnnlib/cGAN_imgSz64_nCh_1_prop_0_input_eps_0.010_output_eps_0.015.vnnlib";
"onnx/cGAN_imgSz64_nCh_1.onnx",	"vnnlib/cGAN_imgSz64_nCh_1_prop_1_input_eps_0.010_output_eps_0.015.vnnlib";
"onnx/cGAN_imgSz64_nCh_1.onnx",	"vnnlib/cGAN_imgSz64_nCh_1_prop_2_input_eps_0.010_output_eps_0.015.vnnlib";
"onnx/cGAN_imgSz64_nCh_1.onnx",	"vnnlib/cGAN_imgSz64_nCh_1_prop_3_input_eps_0.005_output_eps_0.010.vnnlib";
"onnx/cGAN_imgSz64_nCh_3.onnx",	"vnnlib/cGAN_imgSz64_nCh_3_prop_0_input_eps_0.010_output_eps_0.015.vnnlib";
"onnx/cGAN_imgSz64_nCh_3.onnx",	"vnnlib/cGAN_imgSz64_nCh_3_prop_1_input_eps_0.005_output_eps_0.010.vnnlib";
"onnx/cGAN_imgSz64_nCh_3.onnx",	"vnnlib/cGAN_imgSz64_nCh_3_prop_2_input_eps_0.005_output_eps_0.010.vnnlib";
"onnx/cGAN_imgSz64_nCh_3.onnx",	"vnnlib/cGAN_imgSz64_nCh_3_prop_3_input_eps_0.010_output_eps_0.015.vnnlib";
"onnx/cGAN_imgSz32_nCh_3_nonlinear_activations.onnx",	"vnnlib/cGAN_imgSz32_nCh_3_nonlinear_activations_prop_0_input_eps_0.015_output_eps_0.020.vnnlib";
"onnx/cGAN_imgSz32_nCh_1_transposedConvPadding_1.onnx",	"vnnlib/cGAN_imgSz32_nCh_1_transposedConvPadding_1_prop_0_input_eps_0.010_output_eps_0.015.vnnlib";
"onnx/cGAN_imgSz32_nCh_3_upsample.onnx",	"vnnlib/cGAN_imgSz32_nCh_3_upsample_prop_0_input_eps_0.010_output_eps_0.015.vnnlib";
"onnx/cGAN_imgSz32_nCh_3_small_transformer.onnx",	"vnnlib/cGAN_imgSz32_nCh_3_small_transformer_prop_0_input_eps_0.010_output_eps_0.015.vnnlib";
"onnx/cGAN_imgSz32_nCh_3_small_transformer.onnx",	"vnnlib/cGAN_imgSz32_nCh_3_small_transformer_prop_1_input_eps_0.010_output_eps_0.015.vnnlib"    
};
onnx = cgan_instances(i,1);
onnx = onnx{1};
vnnlib = cgan_instances(i,2);
vnnlib = vnnlib{1};
run_vnncomp2025("cgan",onnx,vnnlib,"cgan_results_" + string(i)+".txt");

rmpath(genpath(dir))

end