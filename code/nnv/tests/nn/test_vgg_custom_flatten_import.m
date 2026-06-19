function tests = test_vgg_custom_flatten_import
% Regression test for matlab2nnv support of importer-generated CUSTOM flatten layers.
%
% The MATLAB ONNX importer emits per-model custom Formattable flatten layers (e.g.
% vgg16_7.FlattenLayer1000) for the no-op Flatten ops ONNX inserts before each Gemm. matlab2nnv
% had no mapping for these -> "Unsupported Class of Layer", which blocked the official VNN-COMP
% vgg16-7.onnx. matlab2nnv now maps custom (non-'nnet.') flatten classes to NNV's FlattenLayer
% (identity on already-flat input). This test imports the official VGG16, converts it, and
% CROSS-VALIDATES that NN.evaluate matches the imported dlnetwork's predict (translation
% correctness) -- a custom flatten that was NOT a true no-op would change element order and fail.
%
% Requires vgg16-7.onnx (553 MB, not in the repo): set env NNV_VGG_ONNX, or it is auto-detected at
% the lambda-box benchmark path. The test SKIPS (assumeTrue) when the onnx is unavailable, so it is
% a no-op in environments without the model.
%
% Run: runtests('test_vgg_custom_flatten_import')
    tests = functiontests(localfunctions);
end

function onnx = i_locate_vgg_onnx()
    onnx = getenv('NNV_VGG_ONNX');
    if ~isempty(onnx) && isfile(onnx), return; end
    cands = {
        '/home/johnsott/vsc_nnv-2026/vnncomp2026_benchmarks/benchmarks/vggnet16_2022/2.0/onnx/vgg16-7.onnx'
        '/home/johnsott/vsc_nnv-2026/vnncomp2026_benchmarks/benchmarks/vggnet16_2022/1.0/onnx/vgg16-7.onnx'
        };
    onnx = '';
    for k = 1:numel(cands)
        if isfile(cands{k}), onnx = cands{k}; return; end
    end
end

function test_matlab2nnv_imports_custom_flatten(tc)
    onnx = i_locate_vgg_onnx();
    tc.assumeTrue(~isempty(onnx) && isfile(onnx), ...
        'vgg16-7.onnx not available -> skip (set NNV_VGG_ONNX to the official onnx path)');

    net = importNetworkFromONNX(onnx, "InputDataFormats","BCSS", "OutputDataFormats","BC");

    % matlab2nnv must NOT throw "Unsupported Class of Layer" on the custom flatten layers.
    nnvnet = matlab2nnv(net);
    tc.verifyClass(nnvnet, 'NN');

    % TRANSLATION CORRECTNESS: NN.evaluate must match the imported dlnetwork's predict. A custom
    % flatten mapped with the wrong element order would diverge here.
    x = rand([224 224 3], 'single');
    yref = predict(net, dlarray(x, 'SSCB'));
    yref = double(gather(extractdata(yref))); yref = yref(:);
    ynnv = double(nnvnet.evaluate(double(x))); ynnv = ynnv(:);

    tc.verifyEqual(numel(ynnv), numel(yref), 'output size mismatch');
    tc.verifyLessThan(max(abs(ynnv - yref)), 1e-2, ...
        'NN.evaluate must match dlnetwork.predict (custom-flatten translation correctness)');
    [~, amr] = max(yref); [~, amn] = max(ynnv);
    tc.verifyEqual(amn, amr, 'argmax must agree between NN.evaluate and dlnetwork.predict');
end
