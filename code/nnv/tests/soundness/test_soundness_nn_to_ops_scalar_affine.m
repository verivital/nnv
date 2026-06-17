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
nIn = 6;   % lsnc relu_quadrotor2d_state: flat [6] feature input

% (1) WITHOUT the input dim (the old call signature): a scalar ElementwiseAffine on the flat input
% is still refused with elemAffineFlatSize (the fix is opt-in via inputDim; no behavior change).
oldId = '';
try, nn_to_ops(net, 'colmajor'); catch e, oldId = e.identifier; end
assert(strcmp(oldId, 'nn_to_ops:elemAffineFlatSize'), ...
    'without inputDim, the scalar ElementwiseAffine should still refuse (elemAffineFlatSize), got "%s"', oldId);

% (2) WITH the input dim: the scalar ElementwiseAffine is RESOLVED; nn_to_ops advances past it and
% only refuses at a genuinely-later unsupported op (lsnc has Concat + bilinear ElementwiseProduct).
newId = ''; newMsg = '';
try, nn_to_ops(net, 'colmajor', nIn); catch e, newId = e.identifier; newMsg = e.message; end
assert(~isempty(newId), 'lsnc is not fully supported yet (Concat/bilinear) -> expected a later refusal');
assert(~strcmp(newId, 'nn_to_ops:elemAffineFlatSize'), ...
    'with inputDim, the scalar ElementwiseAffine must be RESOLVED, not refused; got elemAffineFlatSize again');
assert(contains(newMsg,'Concat') || contains(newMsg,'multi-input') || contains(newMsg,'Product') || contains(newMsg,'unsupported'), ...
    'expected advance to a later op (Concat / ElementwiseProduct), got: %s', newMsg);

fprintf('test_soundness_nn_to_ops_scalar_affine: passed (scalar affine resolved with inputDim; advances to "%s")\n', newId);
