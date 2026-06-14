% test_soundness_gpu_bab_maxpool
% Soundness tests for the GPU-BaB engine's MaxPooling2D support. MaxPool is the only
% NONLINEAR op besides ReLU: monotone -> exact IBP; CROWN uses a sound per-window
% relaxation (lower bound = the window's argmax-LOWER input; upper bound = that input if it
% DECIDES the max, else the window's max-upper constant). Validated on a conv + ReLU +
% non-overlapping maxpool + conv + ReLU + OVERLAPPING maxpool + FC net.
% Each %% section runs independently under runtests; net is kept small for speed.

rng(3);
layers = [imageInputLayer([8 8 2],'Normalization','none','Name','in')
          convolution2dLayer(3,4,'Padding',1,'Name','c1'); reluLayer('Name','r1')
          maxPooling2dLayer(2,'Stride',2,'Name','mp1')             % 8 -> 4 (non-overlap)
          convolution2dLayer(3,3,'Padding',1,'Name','c2'); reluLayer('Name','r2')
          maxPooling2dLayer(3,'Stride',1,'Name','mp2')             % 4 -> 2 (OVERLAPPING)
          fullyConnectedLayer(5,'Name','fc')];
nnvnet = matlab2nnv(dlnetwork(layers));
ops    = nn_to_ops(nnvnet);
Cspec  = [eye(4) -ones(4,1)];
xc     = rand(8,8,2);
lb     = xc(:) - 0.04;  ub = xc(:) + 0.04;
Xs     = lb + (ub - lb) .* rand(128, 2000);
tmp = inf(4,1);
for s = 1:size(Xs,2)
    ys = nnvnet.evaluate(reshape(Xs(:,s), [8 8 2]));
    tmp = min(tmp, Cspec * ys(:));
end
trueMin = tmp;

%% Test 1: maxpool ops extracted; degenerate IBP == NNV evaluate
xt = rand(8,8,2); yn = nnvnet.evaluate(xt); yn = yn(:);
yo = gpu_bab_ibp(ops, xt(:), xt(:), 'double');
assert(max(abs(yo(:) - yn)) < 1e-4, 'degenerate maxpool IBP must match NNV evaluate');

%% Test 2: both maxpool ops present (incl. the overlapping one)
types = cellfun(@(o) string(o.type), ops);
assert(sum(types == "maxpool") == 2, 'expected 2 maxpool ops (one overlapping stride<pool)');

%% Test 3: IBP box is sound (contains every sampled output)
[olb, oub] = gpu_bab_ibp(ops, lb, ub, 'double');
ok = true;
for s = 1:size(Xs,2)
    ys = nnvnet.evaluate(reshape(Xs(:,s), [8 8 2])); ys = ys(:);
    if any(ys < olb - 1e-6) || any(ys > oub + 1e-6), ok = false; break; end
end
assert(ok, 'maxpool IBP box must contain all sampled outputs');

%% Test 4: degenerate CROWN is exact (every window decided at a point)
xt = rand(8,8,2); md = gpu_bab_crown_tight(ops, xt(:), xt(:), Cspec, 'double');
yn = nnvnet.evaluate(xt);
assert(max(abs(md - Cspec * yn(:))) < 1e-4, 'degenerate maxpool CROWN must equal C*evaluate');

%% Test 5: CROWN-tight is sound and no looser than IBP
mt = gpu_bab_crown_tight(ops, lb, ub, Cspec, 'double');
assert(all(mt <= trueMin + 1e-4), 'maxpool CROWN-tight margin must be sound (<= true min)');
[ilb, iub] = gpu_bab_ibp(ops, lb, ub, 'double');
Cp = max(Cspec, 0); Cn = min(Cspec, 0);
assert(all(mt >= Cp*ilb + Cn*iub - 1e-6), 'maxpool CROWN-tight must be no looser than IBP');

%% Summary
disp('test_soundness_gpu_bab_maxpool: all sections passed');
