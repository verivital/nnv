function tests = test_matlab2nnv_reshape_vars
%TEST_MATLAB2NNV_RESHAPE_VARS  Regression test for the R2026a fused-Reshape import gap.
%
%   importNetworkFromONNX (R2026a) emits custom layers that carry their ONNX
%   initializers in a `.Vars` property (no `.ONNXParams`). matlab2nnv's fused
%   double-reshape no-op branch only recognized the old ONNXParams form, so a
%   `Reshape_To_ReshapeLayer*` custom layer fell through to ReshapeLayer.parse
%   and died on `Unrecognized ... 'ONNXParams'` -- which silently turned ~50
%   provable dist_shift_2023 UNSATs into `unknown` in the VNN-COMP sweep (the
%   runner's is_nnvnet_valid skip). The fix accepts the .Vars form as a no-op
%   ONLY when it is a 2-shape chain whose final target is flat ([-1 N]), i.e.
%   reshape-to-image-then-flatten-back = identity on flat data.
%
%   The integration test needs the dist_shift_2023 model from the local
%   vnncomp2026_benchmarks clone (or NNV_VNNCOMP2026_BENCHMARKS env var) and is
%   skipped via assumeTrue when unavailable, so CI stays green without it.
%
%   To run: results = runtests('test_matlab2nnv_reshape_vars')
    tests = functiontests(localfunctions);
end

function test_mnist_concat_converts_and_is_faithful(tc)
    onnx = locate_mnist_concat(tc);
    w = warning('off', 'all');
    cleaner = onCleanup(@() warning(w));
    net = importNetworkFromONNX(onnx, "InputDataFormats","BC", "OutputDataFormats","BC");

    % the import must actually contain the R2026a .Vars-form fused reshape this
    % test pins -- if MathWorks changes the import shape again, surface it here
    idx = find(arrayfun(@(L) contains(class(L), 'Reshape_To_ReshapeLayer'), net.Layers), 1);
    assumeTrue(tc, ~isempty(idx), 'import no longer emits a Reshape_To_ReshapeLayer custom layer');
    assumeTrue(tc, isprop(net.Layers(idx), 'Vars'), 'custom layer no longer uses the .Vars form');

    % pre-fix this errored: Unrecognized ... 'ONNXParams' for class '...Reshape_To_ReshapeLayer1000'
    nnvnet = matlab2nnv(net);
    verifyEqual(tc, numel(nnvnet.Layers), numel(net.Layers), ...
        'every MATLAB layer maps to one NNV layer (fused reshape becomes a placeholder no-op)');

    % conversion fidelity: NNV evaluate must match dlnetwork predict (this is
    % what makes downstream approx-star UNSAT proofs trustworthy)
    rng(7);
    % dlarray label note: "InputDataFormats","BC" describes the ONNX tensor layout
    % (batch x channel); the imported dlnetwork has a featureInputLayer(792), and a
    % MATLAB column vector x [792x1] is labeled per-DIMENSION as 'CB' (C=792
    % features, B=1 batch). This orientation is empirically pinned: the same
    % forward pass matches onnxruntime on the original ONNX to ~2e-6 relative
    % (2026-06-12 xval) -- a transposed feed could not reproduce that.
    maxAbs = 0; maxMag = 0;
    for k = 1:5
        x = randn(792, 1, 'single');
        y_nnv = nnvnet.evaluate(double(x));
        y_dl  = extractdata(predict(net, dlarray(x, 'CB')));
        maxAbs = max(maxAbs, max(abs(double(y_nnv(:)) - double(y_dl(:)))));
        maxMag = max(maxMag, max(abs(double(y_dl(:)))));
    end
    relErr = maxAbs / max(maxMag, 1);
    verifyLessThan(tc, relErr, 1e-3, ...
        sprintf('NNV forward pass must match dlnetwork (rel err %.2e; competition tolerance 1e-3)', relErr));
end

% ---------- helpers ----------

function onnx = locate_mnist_concat(tc)
    root = getenv('NNV_VNNCOMP2026_BENCHMARKS');
    if isempty(root)
        here = fileparts(mfilename('fullpath'));
        root = fullfile(here, '..', '..', '..', '..', '..', 'vnncomp2026_benchmarks');
    end
    if ~isfolder(fullfile(root, 'benchmarks')) && isfolder(root)
        % env var may point directly at the benchmarks dir
        cand = root;
    else
        cand = fullfile(root, 'benchmarks');
    end
    onnx = fullfile(cand, 'dist_shift_2023', '1.0', 'onnx', 'mnist_concat.onnx');
    if ~isfile(onnx) && isfile([onnx '.gz'])
        gunzip([onnx '.gz']);
    end
    assumeTrue(tc, isfile(onnx), ...
        sprintf('dist_shift model not found at %s -- skipping (needs the local benchmarks clone)', onnx));
end
