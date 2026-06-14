% test_soundness_gpu_bab_conv
% Soundness tests for the GPU-BaB engine's CONV layer support (Conv2D + folded
% ImageInputLayer normalization). Validates the dlconv IBP (W+/W- split) and the EXACT
% dltranspconv CROWN adjoint on a net with padding + stride + normalization -- the cases
% that exercise the conv<->matrix reshape, the input-ordering, and the strided adjoint.
% To run: results = runtests('test_soundness_gpu_bab_conv')
%
% NOTE: each %% section runs independently under runtests, so every shared quantity is
% computed here in the shared-variables setup (not in a prior test).

% ---- shared setup (deterministic): conv net (pad+stride+norm), ops, box, MC truth ----
rng(1);
layers = [imageInputLayer([10 10 2],'Normalization','zerocenter','Mean',0.3*ones(10,10,2),'Name','in')
          convolution2dLayer(3,5,'Padding',1,'Stride',1,'Name','c1'); reluLayer('Name','r1')
          convolution2dLayer(3,3,'Padding',1,'Stride',2,'Name','c2'); reluLayer('Name','r2')
          fullyConnectedLayer(4,'Name','fc')];
nnvnet = matlab2nnv(dlnetwork(layers));
ops    = nn_to_ops(nnvnet);
Cspec  = [1 -1 0 0; 1 0 -1 0; 1 0 0 -1];           % robustness spec for class 1
xc     = rand(10,10,2);
lb     = xc(:) - 0.05;  ub = xc(:) + 0.05;
Xs     = lb + (ub - lb) .* rand(200, 3000);        % Monte-Carlo input samples
tmp = inf(3,1);
for s = 1:size(Xs,2)
    ys = nnvnet.evaluate(reshape(Xs(:,s), [10 10 2]));
    tmp = min(tmp, Cspec * ys(:));
end
trueMin = tmp;

%% Test 1: conv op-extraction + input-ordering guard (degenerate IBP == NNV evaluate)
xt = rand(10,10,2);
yn = nnvnet.evaluate(xt); yn = yn(:);
yo = gpu_bab_ibp(ops, xt(:), xt(:), 'double');
assert(max(abs(yo(:) - yn)) < 1e-4, 'degenerate-box conv IBP must match NNV evaluate (input ordering)');

%% Test 2: conv IBP box is sound (contains every sampled output)
[olb, oub] = gpu_bab_ibp(ops, lb, ub, 'double');
ok = true;
for s = 1:size(Xs,2)
    ys = nnvnet.evaluate(reshape(Xs(:,s), [10 10 2])); ys = ys(:);
    if any(ys < olb - 1e-6) || any(ys > oub + 1e-6), ok = false; break; end
end
assert(ok, 'conv IBP box must contain all sampled outputs');

%% Test 3: conv CROWN backward is EXACT at a point (== C*evaluate, the dltranspconv adjoint)
xt = rand(10,10,2);
md = gpu_bab_crown_tight(ops, xt(:), xt(:), Cspec, 'double');
yn = nnvnet.evaluate(xt);
assert(max(abs(md - Cspec * yn(:))) < 1e-4, 'degenerate conv CROWN must equal C*evaluate (exact adjoint)');

%% Test 4: conv CROWN-tight is sound and no looser than IBP
mt = gpu_bab_crown_tight(ops, lb, ub, Cspec, 'double');
assert(all(mt <= trueMin + 1e-4), 'conv CROWN-tight margin must be sound (<= true min)');
[ilb, iub] = gpu_bab_ibp(ops, lb, ub, 'double');
Cp = max(Cspec, 0); Cn = min(Cspec, 0);
assert(all(mt >= Cp*ilb + Cn*iub - 1e-6), 'conv CROWN-tight must be no looser than IBP');

%% Summary
disp('test_soundness_gpu_bab_conv: all sections passed');
