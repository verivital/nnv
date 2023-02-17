%% Examples to load nn4sys neural networks

% path
vnnFolder = "/home/manzand/Documents/MATLAB/vnncomp2022_benchmarks/benchmarks/";
nn4sysFolder = "nn4sys/onnx/";

% lindex
net = onnx2nnv(vnnFolder + nn4sysFolder + "lindex.onnx");

% lindex_deep
net1 = onnx2nnv(vnnFolder + nn4sysFolder + "lindex_deep.onnx");

% The remaining networks contain slice layers so we cannot support them
