% test_soundness_gpu_bab_bn_pool
% Soundness tests for the GPU-BaB engine's BatchNorm + AveragePooling layer support.
% BatchNorm (inference) folds to a per-channel affine; average/global-average pooling is
% a per-channel uniform (monotone) linear map. Both are LINEAR -> exact CROWN, no
% relaxation. Validated on a conv + BN + ReLU + avgpool + conv + BN + ReLU + GLOBAL-avgpool
% + FC net (a resnet-style feedforward backbone), with non-trivial injected BN statistics.
% To run: results = runtests('test_soundness_gpu_bab_bn_pool')
%
% NOTE: each %% section runs independently under runtests, so every shared quantity is
% computed here in the shared-variables setup.

% ---- shared setup: build the net, inject non-trivial BN stats, extract ops, MC truth ----
rng(2);
layers = [imageInputLayer([16 16 3],'Normalization','zscore', ...
              'Mean',0.4*ones(16,16,3),'StandardDeviation',0.2*ones(16,16,3),'Name','in')
          convolution2dLayer(3,6,'Padding',1,'Stride',1,'Name','c1')
          batchNormalizationLayer('Name','bn1'); reluLayer('Name','r1')
          averagePooling2dLayer(2,'Stride',2,'Name','ap1')             % 16 -> 8
          convolution2dLayer(3,4,'Padding',1,'Stride',2,'Name','c2')   % 8 -> 4
          batchNormalizationLayer('Name','bn2'); reluLayer('Name','r2')
          globalAveragePooling2dLayer('Name','gap')                    % 4x4 -> 1x1
          fullyConnectedLayer(5,'Name','fc')];
net = dlnetwork(layers);
S = net.State;                              % inject non-trivial BN running stats
for i = 1:height(S)
    if strcmp(S.Parameter{i},'TrainedMean'),     S.Value{i} = dlarray(0.3*randn(size(S.Value{i}),'single')); end
    if strcmp(S.Parameter{i},'TrainedVariance'), S.Value{i} = dlarray(0.5+0.5*rand(size(S.Value{i}),'single')); end
end
net.State = S;
Lr = net.Learnables;                        % inject non-trivial BN scale/offset
for i = 1:height(Lr)
    if contains(Lr.Layer{i},'bn')
        if strcmp(Lr.Parameter{i},'Scale'),  Lr.Value{i} = dlarray(0.5+0.5*rand(size(Lr.Value{i}),'single')); end
        if strcmp(Lr.Parameter{i},'Offset'), Lr.Value{i} = dlarray(0.2*randn(size(Lr.Value{i}),'single')); end
    end
end
net.Learnables = Lr;

nnvnet = matlab2nnv(net);
ops    = nn_to_ops(nnvnet);
Cspec  = [eye(4) -ones(4,1)];               % class-5 robustness margins out(c)-out(5)
xc     = rand(16,16,3);
lb     = xc(:) - 0.03;  ub = xc(:) + 0.03;
Xs     = lb + (ub - lb) .* rand(768, 3000);
tmp = inf(4,1);
for s = 1:size(Xs,2)
    ys = nnvnet.evaluate(reshape(Xs(:,s), [16 16 3]));
    tmp = min(tmp, Cspec * ys(:));
end
trueMin = tmp;

%% Test 1: BN folds + avgpool extract; degenerate IBP == NNV evaluate (input ordering)
xt = rand(16,16,3);
yn = nnvnet.evaluate(xt); yn = yn(:);
yo = gpu_bab_ibp(ops, xt(:), xt(:), 'double');
assert(max(abs(yo(:) - yn)) < 1e-4, 'degenerate-box IBP must match NNV evaluate (BN fold + avgpool + ordering)');

%% Test 2: BN + avgpool ops appear in the decomposition (no silent skip)
types = cellfun(@(o) string(o.type), ops);
assert(sum(types=="normaffine") >= 3, 'expected input-norm + 2 BatchNorm folds as normaffine ops');
assert(sum(types=="avgpool") == 2, 'expected the averagePooling + globalAveragePooling ops');

%% Test 3: IBP box is sound (contains every sampled output)
[olb, oub] = gpu_bab_ibp(ops, lb, ub, 'double');
ok = true;
for s = 1:size(Xs,2)
    ys = nnvnet.evaluate(reshape(Xs(:,s), [16 16 3])); ys = ys(:);
    if any(ys < olb - 1e-6) || any(ys > oub + 1e-6), ok = false; break; end
end
assert(ok, 'IBP box must contain all sampled outputs (BN + avgpool)');

%% Test 4: CROWN backward is EXACT at a point (BN backward + avgpool repelem adjoint)
xt = rand(16,16,3);
md = gpu_bab_crown_tight(ops, xt(:), xt(:), Cspec, 'double');
yn = nnvnet.evaluate(xt);
assert(max(abs(md - Cspec * yn(:))) < 1e-4, 'degenerate CROWN must equal C*evaluate (exact BN + avgpool adjoints)');

%% Test 5: CROWN-tight is sound and no looser than IBP
mt = gpu_bab_crown_tight(ops, lb, ub, Cspec, 'double');
assert(all(mt <= trueMin + 1e-4), 'CROWN-tight margin must be sound (<= true min)');
[ilb, iub] = gpu_bab_ibp(ops, lb, ub, 'double');
Cp = max(Cspec, 0); Cn = min(Cspec, 0);
assert(all(mt >= Cp*ilb + Cn*iub - 1e-6), 'CROWN-tight must be no looser than IBP');

%% Summary
disp('test_soundness_gpu_bab_bn_pool: all sections passed');
