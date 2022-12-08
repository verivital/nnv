%% Examples to load cifar2020 benchmarks

% path
vnnFolder = "/home/manzand/Documents/MATLAB/vnncomp2022_benchmarks/benchmarks/";
benchFolder = "cifar2020/onnx/";

% cifar10_2_255_simplified
net = onnx2nnv(vnnFolder+benchFolder+"cifar10_2_255_simplified.onnx");

% cifar10_8_255_simplified
net2 = onnx2nnv(vnnFolder+benchFolder+"cifar10_8_255_simplified.onnx");

% cifar10_2_255
% net1 = onnx2nnv(vnnFolder+benchFolder+"cifar10_2_255.onnx", loadOptions);
% % custom layer, not good need to further look into this problem, but will
%   be very diffucult to come up with a systematic fix

% cifar10_8_255
% net3 = onnx2nnv(vnnFolder+benchFolder+"cifar10_8_255.onnx"); % same as the previous one

%convBigRELU_PGD
net4 = onnx2nnv(vnnFolder+benchFolder+"convBigRELU__PGD.onnx");