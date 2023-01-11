%% Examples to load cifar2020 benchmarks

% path
vnnFolder = "/home/manzand/Documents/MATLAB/vnncomp2022_benchmarks/benchmarks/";
benchFolder = "cifar_biasfield/onnx/";

% base
net = onnx2nnv(vnnFolder+benchFolder+"cifar_base.onnx");

% field 0
loadOptions = struct;
loadoptions.InputDataFormat = "BC";
net2 = onnx2nnv(vnnFolder+benchFolder+"cifar_bias_field_0.onnx", loadOptions);
%  this almost works, we need to add some parsing to the reshape layers in
%  the case that they have a single Nonlearnable parameter like these ones,
%  that specify the reshaping of the input
%  net.Layers(4, 1).ONNXParams.Nonlearnables.x20

% field 1
net1 = onnx2nnv(vnnFolder+benchFolder+"cifar_bias_field_1.onnx");

% field 2
net3 = onnx2nnv(vnnFolder+benchFolder+"cifar_bias_field_2.onnx"); 

% field 10
net4 = onnx2nnv(vnnFolder+benchFolder+"cifar_bias_field_10.onnx");