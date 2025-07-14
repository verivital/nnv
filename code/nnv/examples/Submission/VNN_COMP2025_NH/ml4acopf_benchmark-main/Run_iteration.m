function Run_iteration(i)


%% ml4acopf

disp("Running ml4acopf..")
dir =  pwd;
addpath(genpath(dir));

parentDir = fileparts(dir);
src_dir = parentDir + "\src";
addpath(genpath(src_dir))

ml4acopf_instances =  {
    "onnx/14_ieee_ml4acopf.onnx", "vnnlib/14_ieee_prop1.vnnlib";
    "onnx/14_ieee_ml4acopf.onnx", "vnnlib/14_ieee_prop3.vnnlib";
    "onnx/14_ieee_ml4acopf.onnx", "vnnlib/14_ieee_prop11.vnnlib";
    "onnx/14_ieee_ml4acopf.onnx", "vnnlib/14_ieee_prop7.vnnlib";
    "onnx/14_ieee_ml4acopf.onnx", "vnnlib/14_ieee_prop5.vnnlib";
    "onnx/14_ieee_ml4acopf.onnx", "vnnlib/14_ieee_prop13.vnnlib";
    "onnx/14_ieee_ml4acopf.onnx", "vnnlib/14_ieee_prop9.vnnlib";
    "onnx/14_ieee_ml4acopf.onnx", "vnnlib/14_ieee_prop2.vnnlib";
    "onnx/14_ieee_ml4acopf.onnx", "vnnlib/14_ieee_prop14.vnnlib";
    "onnx/14_ieee_ml4acopf.onnx", "vnnlib/14_ieee_prop10.vnnlib";
    "onnx/14_ieee_ml4acopf.onnx", "vnnlib/14_ieee_prop6.vnnlib";
    "onnx/14_ieee_ml4acopf.onnx", "vnnlib/14_ieee_prop4.vnnlib";
    "onnx/14_ieee_ml4acopf.onnx", "vnnlib/14_ieee_prop12.vnnlib";
    "onnx/14_ieee_ml4acopf.onnx", "vnnlib/14_ieee_prop8.vnnlib";
    "onnx/300_ieee_ml4acopf.onnx", "vnnlib/300_ieee_prop4.vnnlib";
    "onnx/300_ieee_ml4acopf.onnx", "vnnlib/300_ieee_prop2.vnnlib";
    "onnx/300_ieee_ml4acopf.onnx", "vnnlib/300_ieee_prop105.vnnlib";
    "onnx/300_ieee_ml4acopf.onnx", "vnnlib/300_ieee_prop3.vnnlib";
    "onnx/118_ieee_ml4acopf.onnx", "vnnlib/118_ieee_prop6.vnnlib";
    "onnx/118_ieee_ml4acopf.onnx", "vnnlib/118_ieee_prop4.vnnlib";
    "onnx/118_ieee_ml4acopf.onnx", "vnnlib/118_ieee_prop2.vnnlib";
    "onnx/118_ieee_ml4acopf.onnx", "vnnlib/118_ieee_prop5.vnnlib";
    "onnx/118_ieee_ml4acopf.onnx", "vnnlib/118_ieee_prop3.vnnlib";
    };

onnx = ml4acopf_instances(i,1);
onnx = onnx{1};
vnnlib = ml4acopf_instances(i,2);
vnnlib = vnnlib{1};
run_vnncomp2025("ml4acopf",onnx,vnnlib,"ML4aCOPF_results_" + string(i)+".txt");

rmpath(genpath(dir))
end
