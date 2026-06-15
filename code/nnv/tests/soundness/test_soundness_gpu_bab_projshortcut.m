% test_soundness_gpu_bab_projshortcut
% Projection-shortcut (conv-on-skip) DAG soundness for the GPU-BaB FULL-DAG engine. Unlike a
% direct identity skip (test_soundness_gpu_bab_residual), the skip branch here carries a 1x1
% conv whose input is the BLOCK INPUT (the stem r0), several ops BEFORE the add -- NOT the
% immediately-preceding op. master's chain extractor would feed that conv the wrong tensor
% (the previous op's output); the full-DAG per-op `src` routing resolves it from Connections.
% The degenerate-CROWN == C*evaluate check pins the DAG backward routing exactly. This is the
% synthetic analogue of the VNN-COMP cifar100 CIFAR100_resnet_medium block.
% Each %% section runs independently; net kept small.

rng(3);
lg = layerGraph([
  imageInputLayer([8 8 3],'Normalization','none','Name','in')
  convolution2dLayer(3,5,'Padding',1,'Name','c0'); reluLayer('Name','r0')
  convolution2dLayer(3,5,'Padding',1,'Name','c1'); reluLayer('Name','r1')
  convolution2dLayer(3,5,'Padding',1,'Name','c2')
  convolution2dLayer(1,5,'Name','cskip')               % skip-branch projection (1x1 conv)
  additionLayer(2,'Name','add'); reluLayer('Name','r2')
  globalAveragePooling2dLayer('Name','gap')
  fullyConnectedLayer(4,'Name','fc')]);
% Rewire so cskip is a true projection shortcut: its input is the STEM r0 (the block input),
% NOT c2 (the op placed before it in the array); the main branch c2 feeds add/in1.
lg = disconnectLayers(lg, 'c2', 'cskip');
lg = disconnectLayers(lg, 'cskip', 'add/in1');
lg = connectLayers(lg, 'r0', 'cskip');         % skip from the block input (src != previous op)
lg = connectLayers(lg, 'c2', 'add/in1');       % main branch -> add first input
lg = connectLayers(lg, 'cskip', 'add/in2');    % skip branch -> add second input
nnvnet = matlab2nnv(dlnetwork(lg));
ops    = nn_to_ops(nnvnet);
Cspec  = [eye(3) -ones(3,1)];
xc     = rand(8,8,3); lb = xc(:)-0.03; ub = xc(:)+0.03;
Xs     = lb + (ub-lb).*rand(192, 2000);
tmp = inf(3,1);
for s = 1:size(Xs,2), ys = nnvnet.evaluate(reshape(Xs(:,s),[8 8 3])); tmp = min(tmp, Cspec*ys(:)); end
trueMin = tmp;

%% Test 1: a conv-on-skip op routes to a NON-previous src; degenerate IBP == NNV evaluate
types = cellfun(@(o) string(o.type), ops);
assert(sum(types=="add")==1, 'expected one add op');
assert(any(arrayfun(@(k) isfield(ops{k},'src') && ops{k}.src ~= k-1, 1:numel(ops))), ...
    'expected at least one op whose src ~= k-1 (the conv-on-skip projection branch)');
xt = rand(8,8,3); yn = nnvnet.evaluate(xt); yn = yn(:);
yo = gpu_bab_ibp(ops, xt(:), xt(:), 'double');
assert(max(abs(yo(:)-yn)) < 1e-4, 'degenerate projection-shortcut IBP must match NNV evaluate');

%% Test 2: IBP box is sound (contains every sampled output)
[olb,oub] = gpu_bab_ibp(ops, lb, ub, 'double'); ok = true;
for s = 1:size(Xs,2)
    ys = nnvnet.evaluate(reshape(Xs(:,s),[8 8 3])); ys = ys(:);
    if any(ys < olb-1e-6) || any(ys > oub+1e-6), ok = false; break; end
end
assert(ok, 'projection-shortcut IBP box must contain all MC outputs');

%% Test 3: degenerate CROWN is exact (the DAG skipA backward through the conv-on-skip)
xt = rand(8,8,3); md = gpu_bab_crown_tight(ops, xt(:), xt(:), Cspec, 'double'); yn = nnvnet.evaluate(xt);
assert(max(abs(md - Cspec*yn(:))) < 1e-4, 'degenerate projection-shortcut CROWN must equal C*evaluate (DAG skipA exact)');

%% Test 4: CROWN-tight is sound and no looser than IBP
mt = gpu_bab_crown_tight(ops, lb, ub, Cspec, 'double');
assert(all(mt <= trueMin + 1e-4), 'projection-shortcut CROWN-tight must be sound (<= true min)');
[ilb,iub] = gpu_bab_ibp(ops, lb, ub, 'double'); Cp = max(Cspec,0); Cn = min(Cspec,0);
assert(all(mt >= Cp*ilb + Cn*iub - 1e-6), 'projection-shortcut CROWN-tight must be no looser than IBP');

%% Summary
disp('test_soundness_gpu_bab_projshortcut: all sections passed');
