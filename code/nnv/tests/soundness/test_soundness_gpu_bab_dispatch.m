% test_soundness_gpu_bab_dispatch
% Soundness tests for gpu_bab_try_verify -- the additive GPU-BaB pre-check. The cardinal
% property: it must NEVER return 'robust' when the box contains a counterexample (a -150),
% and any 'unsafe' must be a real witness. It must 'skip' unsupported architectures.
% Small nets / modest budgets for speed.

rng(11);
mklayers = @() [imageInputLayer([6 6 1],'Normalization','none','Name','in')
                convolution2dLayer(3,4,'Padding',1,'Name','c1'); reluLayer
                maxPooling2dLayer(2,'Stride',2,'Name','mp')
                fullyConnectedLayer(8,'Name','f1'); reluLayer
                fullyConnectedLayer(3,'Name','fc')];
net  = matlab2nnv(dlnetwork(mklayers()));
opts = struct('precision','double','maxNodes',150,'nSample',24,'cexEvery',10);

%% Test 1: a tiny box that is robust -> 'robust' or 'unknown'; if robust, NO MC counterexample
x = rand(6,6,1); y = net.evaluate(x); [~,tgt] = max(y(:));
lb = max(x(:)-0.015,0); ub = min(x(:)+0.015,1);
[v, inf] = gpu_bab_try_verify(net, lb, ub, tgt, opts);
assert(ismember(string(v), ["robust","unknown"]), 'tiny robust box must not be unsafe/skip');
assert(inf.guardErr < 1e-4, 'orientation guard must pass on a correctly-extracted net');
if strcmp(v,'robust')
    bad = false;
    for s=1:4000, xs=lb+(ub-lb).*rand(36,1); ys=net.evaluate(reshape(xs,[6 6 1])); [~,p]=max(ys(:)); if p~=tgt, bad=true; break; end, end
    assert(~bad, 'CARDINAL: returned robust but MC found a counterexample');
end

%% Test 2: a box that CONTAINS a counterexample -> must NOT be 'robust'
x = rand(6,6,1); y = net.evaluate(x); [~,tgt] = max(y(:));
xp = x; found = false;
for k=1:400
    cand = min(max(x + 0.7*randn(6,6,1),0),1);
    yp = net.evaluate(cand); [~,p] = max(yp(:));
    if p ~= tgt, xp = cand; found = true; break; end
end
assert(found, 'test setup: could not find a flipping point');
lb = min(x(:),xp(:)); ub = max(x(:),xp(:));        % box spans x and the counterexample xp
[v2, ~] = gpu_bab_try_verify(net, lb, ub, tgt, opts);
assert(~strcmp(v2,'robust'), 'CARDINAL: certified robust on a box that contains a counterexample (-150)');

%% Test 3: unsupported architecture (tanh activation) -> 'skip'
ulayers = [featureInputLayer(4,'Name','in'); fullyConnectedLayer(5,'Name','f1'); ...
           tanhLayer('Name','t'); fullyConnectedLayer(3,'Name','out')];
unet = matlab2nnv(dlnetwork(ulayers));
[v3, i3] = gpu_bab_try_verify(unet, -ones(4,1), ones(4,1), 1, opts);
assert(strcmp(v3,'skip'), 'unsupported arch (tanh) must skip, got %s', v3);

%% Test 4: any 'unsafe' verdict carries a witness that truly misclassifies under net.evaluate
% (sweep a few boxes; only assert when 'unsafe' actually fires)
nUnsafe = 0;
for t = 1:6
    xx = rand(6,6,1); yy = net.evaluate(xx); [~,tt] = max(yy(:));
    e = 0.2*t;
    lbi = max(xx(:)-e,0); ubi = min(xx(:)+e,1);
    [vi, ii] = gpu_bab_try_verify(net, lbi, ubi, tt, opts);
    if strcmp(vi,'unsafe')
        nUnsafe = nUnsafe + 1;
        yc = net.evaluate(reshape(ii.cex,[6 6 1])); [~,pc] = max(yc(:));
        assert(pc ~= tt, 'unsafe witness must misclassify under net.evaluate');
    end
end
fprintf('(%d unsafe verdicts, all witness-confirmed)\n', nUnsafe);

%% Test 5: flatten guard distinguishes a WRONG order on a UNIFORM L-inf box (audit -150 fix)
% Worst case: a constant image + uniform-eps box. Collinear probes lb+s*(ub-lb) were BLIND to
% column permutations there ((P-I)*ones=0), so a wrong flatten order could pass and bound the
% wrong function. The non-uniform (low-discrepancy) probes must REJECT the wrong orders.
net2 = matlab2nnv(dlnetwork(mklayers()));
x0 = 0.5*ones(6,6,1); ep = 0.05; lbf = x0(:)-ep; ubf = x0(:)+ep;   % constant image, uniform L-inf
nn = numel(lbf); ii = (1:nn)';
g1 = mod(ii*0.6180339887498949,1); g2 = mod(ii*1.3247179572447460+0.37,1); g3 = mod(ii*0.7548776662466927+0.11,1);
pr = [(lbf+ubf)/2, lbf+g1.*(ubf-lbf), lbf+g2.*(ubf-lbf), lbf+g3.*(ubf-lbf), lbf+(1-g1).*(ubf-lbf)];
ords = {'colmajor','chw_rowmajor','hwc_rowmajor'}; matchOrder = false(1,3);
for oi = 1:3
    ops2 = nn_to_ops(net2, ords{oi}); good = true;
    for pp = 1:size(pr,2)
        cp = pr(:,pp); yo = gpu_bab_ibp(ops2, cp, cp, 'double'); yn = net2.evaluate(reshape(cp,[6 6 1]));
        if max(abs(yo(:)-yn(:))) > 1e-4*max(1,max(abs(yn(:)))), good = false; break; end
    end
    matchOrder(oi) = good;
end
assert(matchOrder(1), 'colmajor (the native flatten) must match at the non-uniform probes');
assert(~matchOrder(2) && ~matchOrder(3), 'CARDINAL: wrong flatten orders must be REJECTED on the uniform-box worst case (else false-robust)');

%% Test 6: a caller-supplied precision='single' must not yield an unsound robust (forced double)
x = rand(6,6,1); y = net.evaluate(x); [~,tgt] = max(y(:));
lb = max(x(:)-0.015,0); ub = min(x(:)+0.015,1);
[vs, ~] = gpu_bab_try_verify(net, lb, ub, tgt, struct('precision','single','maxNodes',120));
if strcmp(vs,'robust')
    bad = false;
    for s=1:3000, xs=lb+(ub-lb).*rand(36,1); ys=net.evaluate(reshape(xs,[6 6 1])); [~,p]=max(ys(:)); if p~=tgt, bad=true; break; end, end
    assert(~bad, 'precision=single robust must still be MC-sound (engine forces double internally)');
end

%% Summary
disp('test_soundness_gpu_bab_dispatch: all sections passed');
