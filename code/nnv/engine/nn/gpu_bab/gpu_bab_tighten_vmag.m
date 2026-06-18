function vmag = gpu_bab_tighten_vmag(vmag, rootBounds, ops, reluIdx)
% GPU_BAB_TIGHTEN_VMAG  Tighten the M3b value-magnitude majorant using crown-tight bounds.
%   vmag = GPU_BAB_TIGHTEN_VMAG(vmag, rootBounds, ops, reluIdx)
%
%   `vmag{k}` is a DOUBLE per-op output value-magnitude majorant used by gpu_bab_crown_spec_dag's
%   per-op `derr` (the sound-FP32 outward widening: derr ~ sum_op gamma(m)*|A|*vmag). Computed from
%   gpu_bab_ibp, it WRAPS and blows up with depth on conv resnets (measured: up to ~500 at deep
%   layers), which inflates `derr` so much the widened margin cannot certify (diagnosed 2026-06-18:
%   idx_8762 binding derr=0.38 vs an un-widened margin the deep BaB does close).
%
%   gpu_bab_crown_tight already computes TIGHT pre-activation bounds (rootBounds.preL/preU) for every
%   ReLU -- far tighter than IBP and SOUND (outward-widened in single). For each ReLU op k:
%     * its INPUT  magnitude is bounded by max(|preL{k}|, |preU{k}|)  (= the op feeding it, vmag{src+1})
%     * its OUTPUT magnitude is bounded by max(0, preU{k})            (post-ReLU is in [0, max(0,preU)])
%   Replacing the IBP value with the elementwise MIN of (IBP, crown-tight) keeps a SOUND over-estimate
%   (both bound |output|; min preserves the bound) while shrinking `derr` toward certifiability.
%
%   SOUNDNESS: min of two sound magnitude majorants is a sound magnitude majorant. The crown-tight
%   bounds are sound lower/upper bounds (outward-widened in single, exact in double), so their
%   abs-max is >= the true |output| over the box. Indices are mapped via ops{k}.src (the op feeding
%   ReLU k; 0 = net input), and a size guard skips any mismatch (keeps the IBP value -> still sound).
%   No-op when vmag/rootBounds absent (the unsound-screen / double paths are untouched).

    if isempty(vmag) || isempty(rootBounds) || ~isfield(rootBounds, 'preL'), return; end
    preL = rootBounds.preL; preU = rootBounds.preU;
    for r = 1:numel(reluIdx)
        k = reluIdx(r);
        if k > numel(preL) || isempty(preL{k}) || ~isfield(ops{k}, 'src'), continue; end
        % relu k's INPUT (= the output of op ops{k}.src; vmag index src+1)
        si = ops{k}.src + 1;
        preMag = double(max(abs(gather(preL{k}(:))), abs(gather(preU{k}(:)))));
        if si >= 1 && si <= numel(vmag) && numel(vmag{si}) == numel(preMag)
            vmag{si} = min(double(vmag{si}(:)), preMag);
        end
        % relu k's OUTPUT (post-ReLU in [0, max(0,preU)]; vmag index k+1)
        outMag = double(max(0, gather(preU{k}(:))));
        if (k + 1) <= numel(vmag) && numel(vmag{k + 1}) == numel(outMag)
            vmag{k + 1} = min(double(vmag{k + 1}(:)), outMag);
        end
    end
end
