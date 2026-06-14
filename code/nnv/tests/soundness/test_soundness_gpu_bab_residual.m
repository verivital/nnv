% test_soundness_gpu_bab_residual
% Residual/DAG (AdditionLayer skip-add) soundness for the GPU-BaB engine. The add is LINEAR
% (out = out[a]+out[b]) -> exact IBP (interval sum) + exact CROWN (the coefficient flows
% UNCHANGED to both inputs via the skipA accumulation). The only risk is the DAG plumbing,
% which the degenerate-CROWN==C*evaluate check pins down. Net: stem -> [conv-relu-conv + skip]
% -> relu -> global-avg-pool -> fc, with the skip from an INTERNAL op (typical resnet block).
% Each %% section runs independently; net kept small.

rng(2);
lg = layerGraph([
  imageInputLayer([8 8 3],'Normalization','none','Name','in')
  convolution2dLayer(3,5,'Padding',1,'Name','c0'); reluLayer('Name','r0')
  convolution2dLayer(3,5,'Padding',1,'Name','c1'); reluLayer('Name','r1')
  convolution2dLayer(3,5,'Padding',1,'Name','c2')
  additionLayer(2,'Name','add'); reluLayer('Name','r2')
  globalAveragePooling2dLayer('Name','gap')
  fullyConnectedLayer(4,'Name','fc')]);
lg = connectLayers(lg, 'r0', 'add/in2');             % skip from r0 (internal op, not input)
nnvnet = matlab2nnv(dlnetwork(lg));
ops    = nn_to_ops(nnvnet);
Cspec  = [eye(3) -ones(3,1)];
xc     = rand(8,8,3); lb = xc(:)-0.03; ub = xc(:)+0.03;
Xs     = lb + (ub-lb).*rand(192, 2000);
tmp = inf(3,1);
for s = 1:size(Xs,2), ys = nnvnet.evaluate(reshape(Xs(:,s),[8 8 3])); tmp = min(tmp, Cspec*ys(:)); end
trueMin = tmp;

%% Test 1: add op extracted (internal skip); degenerate IBP == NNV evaluate
types = cellfun(@(o) string(o.type), ops);
assert(sum(types=="add")==1, 'expected one add op');
xt = rand(8,8,3); yn = nnvnet.evaluate(xt); yn = yn(:);
yo = gpu_bab_ibp(ops, xt(:), xt(:), 'double');
assert(max(abs(yo(:)-yn)) < 1e-4, 'degenerate residual IBP must match NNV evaluate');

%% Test 2: IBP box is sound (contains all MC outputs)
[olb,oub] = gpu_bab_ibp(ops, lb, ub, 'double'); ok = true;
for s = 1:size(Xs,2)
    ys = nnvnet.evaluate(reshape(Xs(:,s),[8 8 3])); ys = ys(:);
    if any(ys < olb-1e-6) || any(ys > oub+1e-6), ok = false; break; end
end
assert(ok, 'residual IBP box must contain all MC outputs');

%% Test 3: degenerate CROWN is exact (the DAG skipA backward)
xt = rand(8,8,3); md = gpu_bab_crown_tight(ops, xt(:), xt(:), Cspec, 'double'); yn = nnvnet.evaluate(xt);
assert(max(abs(md - Cspec*yn(:))) < 1e-4, 'degenerate residual CROWN must equal C*evaluate (skipA exact)');

%% Test 4: CROWN-tight is sound and no looser than IBP
mt = gpu_bab_crown_tight(ops, lb, ub, Cspec, 'double');
assert(all(mt <= trueMin + 1e-4), 'residual CROWN-tight must be sound (<= true min)');
[ilb,iub] = gpu_bab_ibp(ops, lb, ub, 'double'); Cp = max(Cspec,0); Cn = min(Cspec,0);
assert(all(mt >= Cp*ilb + Cn*iub - 1e-6), 'residual CROWN-tight must be no looser than IBP');

%% Test 5: N>2 addition is refused (sound-by-refusal)
lg3 = layerGraph([imageInputLayer([4 4 2],'Normalization','none','Name','in')
  convolution2dLayer(3,2,'Padding',1,'Name','c1')
  additionLayer(3,'Name','add3'); fullyConnectedLayer(2,'Name','fc')]);
lg3 = connectLayers(lg3,'in','add3/in2'); lg3 = connectLayers(lg3,'in','add3/in3');
net3 = matlab2nnv(dlnetwork(lg3));
err = false; try, nn_to_ops(net3); catch, err = true; end
assert(err, 'N=3 addition must be refused for soundness');

%% Summary
disp('test_soundness_gpu_bab_residual: all sections passed');
