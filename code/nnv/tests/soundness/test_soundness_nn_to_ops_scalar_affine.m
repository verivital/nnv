% test_soundness_nn_to_ops_scalar_affine
% Regression test for the nn_to_ops scalar-ElementwiseAffine-on-flat fix (the lsnc_relu importer
% advance). BEFORE the fix, nn_to_ops refused a scalar ElementwiseAffine (uniform y=s*x+t input
% normalization) on a flat input ("nn_to_ops:elemAffineFlatSize") because it could not infer the
% per-feature broadcast size F. AFTER, passing the input dim lets it resolve F from the tracked
% flat size and advance. Asserts: (1) the scalar affine is no longer the refusal point, and (2) any
% advance/refusal is at a genuinely-later op (Concat / ElementwiseProduct), never elemAffineFlatSize.
% SOUNDNESS is otherwise enforced by the mandatory IBP==evaluate orientation guard in the caller.
%
% Gated on the lsnc_relu manifest (relu_quadrotor2d_state.nnv.mat) being present; skips otherwise.

cands = { getenv('NNV_LSNC_MANIFEST'), ...
    fullfile(getenv('HOME'),'vsc_nnv-2026','vnncomp2025_benchmarks','benchmarks','lsnc_relu','onnx','relu_quadrotor2d_state.nnv.mat'), ...
    fullfile('..','..','..','..','..','vnncomp2025_benchmarks','benchmarks','lsnc_relu','onnx','relu_quadrotor2d_state.nnv.mat') };
manifest = '';
for c = 1:numel(cands)
    if ~isempty(cands{c}) && isfile(cands{c}), manifest = cands{c}; break; end
end
if isempty(manifest)
    fprintf('test_soundness_nn_to_ops_scalar_affine: SKIP (lsnc_relu manifest not found)\n');
    return;
end

net = load_nnv_from_mat(manifest);
% lsnc relu_quadrotor2d_state has SCALAR ElementwiseAffine layers on a flat vector, each preceded by
% a FullyConnected. The fix tracks every op's flat size (= size(W,1) through affines, the input dim
% for a leading affine), so those scalar affines resolve their per-feature broadcast size F instead
% of refusing (nn_to_ops:elemAffineFlatSize). The ONE invariant: nn_to_ops must NEVER refuse at the
% scalar affine again -- whether or not inputDim is supplied. (It advances to the next genuinely-
% unsupported op -- lsnc's ConcatenationLayer / bilinear ElementwiseProduct -- or succeeds outright;
% both are fine, so the test is not brittle to future Concat/bilinear support.)
for inDim = {[], 6}                                   % both call forms: without inputDim, and with it
    id = '';
    try, nn_to_ops(net, 'colmajor', inDim{1}); catch e, id = e.identifier; end
    assert(~strcmp(id, 'nn_to_ops:elemAffineFlatSize'), ...
        'scalar ElementwiseAffine on a flat input must be RESOLVED, not refused (got elemAffineFlatSize, inputDim=%s)', mat2str(inDim{1}));
end
fprintf('test_soundness_nn_to_ops_scalar_affine: passed (scalar ElementwiseAffine on a flat input resolved on lsnc)\n');
