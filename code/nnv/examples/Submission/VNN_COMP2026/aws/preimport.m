cd('/home/ubuntu/nnv/code/nnv');
startup_nnv;
od = '/home/ubuntu/vnncomp2025_benchmarks/benchmarks/cifar100_2024/onnx';
xfer = '/tmp/preimport';
if ~isfolder(xfer), mkdir(xfer); end
names = {'CIFAR100_resnet_large', 'CIFAR100_resnet_medium'};
for i = 1:numel(names)
    nm = names{i};
    onnx = fullfile(od, [nm '.onnx']);
    net = importNetworkFromONNX(onnx, "InputDataFormats","BCSS", "OutputDataFormats","BC"); %#ok<NASGU>
    nnvnet = matlab2nnv(net);
    needReshape = 1;                              % cifar100 (confirmed in load_vnncomp_network)
    inShape = nnvnet.Layers{1}.InputSize;
    save(fullfile(xfer, [nm '.mat']), 'nnvnet', 'needReshape', 'inShape', '-v7.3');
    fprintf('SAVED %s inShape=%s\n', nm, mat2str(inShape));
end
disp('PREIMPORT_DONE');
