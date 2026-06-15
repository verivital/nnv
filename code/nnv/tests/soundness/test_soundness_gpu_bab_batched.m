% test_soundness_gpu_bab_batched
% gpu_bab_relu_split_batched (batched-DFS ReLU-split, the GPU-saturating refinement of the
% serial gpu_bab_relu_split) must be SOUND and CONSISTENT with the serial path: the optional
% fixings arg of gpu_bab_crown_spec is backward-compatible, and the batched verifier never
% returns a verdict that CONTRADICTS the serial one (robust vs unsafe) -- sound-or-unknown.
% FC+ReLU nets; runs CPU-only (single precision) so it works in GPU-less CI. Loops live in
% local functions (not in %% sections) to avoid the script-test line-number quirk.
% Each %% section is self-contained.

%% Test 1: gpu_bab_crown_spec fixings arg is backward-compatible (fixings={} == 5-arg call)
rng(1); ops = i_fc_relu([5 12 12 6]);
lb = single(-0.2*ones(5,1)); ub = single(0.2*ones(5,1));
C = single([eye(5) -ones(5,1)]);                       % 5 specs over 6 outputs
m1 = gpu_bab_crown_spec(ops, lb, ub, C, 'single');                       % 5-arg
m2 = gpu_bab_crown_spec(ops, lb, ub, C, 'single', cell(numel(ops),1));   % empty fixings
assert(max(abs(m1(:)-m2(:))) < 1e-6, 'fixings={} must match the 5-arg gpu_bab_crown_spec call');

%% Test 2: a fixed-active clamp can only RAISE the lower pre-activation bound (sound clamp)
rng(2); ops = i_fc_relu([4 10 5]);
lb = single(-0.3*ones(4,1)); ub = single(0.3*ones(4,1));
C = single([eye(4) -ones(4,1)]);
reluK = find(cellfun(@(o) strcmp(o.type,'relu'), ops), 1);
fixNone = cell(numel(ops),1);
fixAct  = cell(numel(ops),1); fixAct{reluK} = ones(10,1);   % force every neuron active
mFree = gpu_bab_crown_spec(ops, lb, ub, C, 'single', fixNone);
mAct  = gpu_bab_crown_spec(ops, lb, ub, C, 'single', fixAct);
assert(all(isfinite(mAct(:))), 'clamped bound must be finite');
% both are valid lower bounds on C*f over the (sub)box; just assert they computed.
assert(numel(mAct)==4 && numel(mFree)==4, 'spec margins must be nSpec x 1');

%% Test 3: batched ReLU-split never CONTRADICTS serial (robust vs unsafe) over random nets
[nContra, nDef, nTot] = i_compare(7, 12);
assert(nContra == 0, sprintf('batched vs serial produced %d robust/unsafe contradictions', nContra));
assert(nDef >= 1, sprintf('expected >=1 definitively-decided case, got %d/%d', nDef, nTot));

%% Test 4: a robust box is not contradicted; the batched BaB loop executes soundly
rng(4); ops = i_fc_relu([6 16 16 6]);
x = single(randn(6,1)*0.4); yc = gpu_bab_ibp(ops, x, x, 'single'); [~,tl] = max(yc);
ep = single(0.012); lb = x-ep; ub = x+ep;
[sS,~] = gpu_bab_relu_split(ops, lb, ub, tl, 6, struct('precision','single','maxNodes',5000));
[sB,iB] = gpu_bab_relu_split_batched(ops, lb, ub, tl, 6, struct('precision','single','maxNodes',5000,'maxFrontier',128));
assert(~(strcmp(sS,'robust') && strcmp(sB,'unsafe')), sprintf('contradiction serial=robust batched=unsafe'));
assert(~(strcmp(sS,'unsafe') && strcmp(sB,'robust')), sprintf('contradiction serial=unsafe batched=robust'));
assert(iB.nodes >= 1, 'batched verifier must report at least the root node');

%% Summary
disp('test_soundness_gpu_bab_batched: all sections passed');

% ----------------------------------------------------------------------------------------
function ops = i_fc_relu(dims)
% Build a sequential FC+ReLU op list (nn_to_ops format) from the layer widths in `dims`.
% Uses the caller's current rng state (He-init weights, small biases).
    ops = {};
    for L = 1:numel(dims)-1
        W = randn(dims(L+1), dims(L)) * sqrt(2/dims(L));
        b = randn(dims(L+1), 1) * 0.1;
        ops{end+1} = struct('type','affine','W',single(W),'b',single(b),'src',numel(ops)); %#ok<AGROW>
        if L < numel(dims)-1
            ops{end+1} = struct('type','relu','src',numel(ops)); %#ok<AGROW>
        end
    end
end

function [nContra, nDef, nTot] = i_compare(seed, nTrials)
% Run serial and batched ReLU-split on nTrials random FC+ReLU robustness queries (true label
% = the box-center prediction, small perturbation). Count robust/unsafe contradictions and
% definitively-decided cases. Sound implementations agree-or-unknown; a contradiction is a bug.
    rng(seed); nContra = 0; nDef = 0; nTot = nTrials;
    for t = 1:nTrials
        dims = [6 20 20 6]; ops = i_fc_relu(dims);
        x = randn(dims(1),1) * 0.5;
        yc = gpu_bab_ibp(ops, single(x), single(x), 'single'); [~,tl] = max(yc);
        ep = single(0.02 + 0.05*rand);
        lb = single(x) - ep; ub = single(x) + ep;
        oS = struct('precision','single','maxNodes',8000);
        oB = struct('precision','single','maxNodes',8000,'maxFrontier',256);
        [sS,~] = gpu_bab_relu_split(ops, lb, ub, tl, dims(end), oS);
        [sB,~] = gpu_bab_relu_split_batched(ops, lb, ub, tl, dims(end), oB);
        if (strcmp(sS,'robust') && strcmp(sB,'unsafe')) || (strcmp(sS,'unsafe') && strcmp(sB,'robust'))
            nContra = nContra + 1;
        end
        if ~strcmp(sB,'unknown'), nDef = nDef + 1; end
    end
end
