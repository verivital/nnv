% test_soundness_gpu_bab_halfspace
% Soundness tests for gpu_bab_halfspace_verify -- the additive general-halfspace (disjunctive)
% GPU-BaB pre-check for the control benchmarks (cersyve / lsnc_relu / linearizenn). The cardinal
% property: it must NEVER return 'robust' when the reachable output ENTERS an unsafe disjunct
% {y : G*y <= g} (a -150). Malformed / missing / mis-dimensioned Hg must 'skip' GRACEFULLY (no
% crash -> the dispatcher falls back to Star). Small FC ReLU control-style net + modest budget.

rng(23);
% control-style net: 4 flat inputs -> 2 outputs (Y_0, Y_1)
layers = [featureInputLayer(4,'Name','in'); fullyConnectedLayer(6,'Name','f1'); ...
          reluLayer('Name','r1'); fullyConnectedLayer(2,'Name','out')];
net = matlab2nnv(dlnetwork(layers));
lb = -ones(4,1); ub = ones(4,1);
opts = struct('maxNodes', 300);                        % bound BaB for test speed

% sample the reachable output range to place disjuncts relative to the reachable set
N = 5000; Y = zeros(2, N);
for k = 1:N, Y(:,k) = net.evaluate(lb + (ub-lb).*rand(4,1)); end
y1 = Y(1,:); med1 = median(y1); min1 = min(y1); max1 = max(y1);

%% Test 1: a disjunct the reachable set clearly AVOIDS -> 'robust' (or 'unknown'); if robust, MC-sound
% unsafe = {y : y_1 <= min1 - 100}: far below every reachable y_1, so provably avoided.
[v1, i1] = gpu_bab_halfspace_verify(net, lb, ub, struct('Hg', HalfSpace([1 0], min1 - 100)), opts);
assert(ismember(string(v1), ["robust","unknown"]), 'avoided disjunct must be robust/unknown, got %s', v1);
assert(i1.guardErr < 1e-4, 'orientation guard must pass on a correctly-extracted FC net');
if strcmp(v1, 'robust')
    bad = false;
    for s = 1:5000, ys = net.evaluate(lb + (ub-lb).*rand(4,1)); if ys(1) <= min1 - 100, bad = true; break; end, end
    assert(~bad, 'CARDINAL: robust but MC found an output inside the unsafe disjunct (-150)');
end

%% Test 2: a disjunct the reachable set ENTERS -> must NOT be 'robust' (-150 guard)
% unsafe = {y : y_1 <= med1}: ~half the reachable outputs are inside, so genuinely unsafe.
[v2, ~] = gpu_bab_halfspace_verify(net, lb, ub, struct('Hg', HalfSpace([1 0], med1)), opts);
assert(~strcmp(v2, 'robust'), 'CARDINAL: certified robust on a disjunct the reachable set enters (-150)');
enters = false;
for s = 1:5000, ys = net.evaluate(lb + (ub-lb).*rand(4,1)); if ys(1) <= med1, enters = true; break; end, end
assert(enters, 'test setup: the unsafe disjunct should be reachable');

%% Test 3: malformed / missing / mis-dimensioned Hg must 'skip' gracefully (the isfield/isprop-safe path)
[va, ~] = gpu_bab_halfspace_verify(net, lb, ub, struct(), opts);                              % no Hg field
assert(strcmp(va, 'skip'), 'missing Hg must skip, got %s', va);
[vb, ~] = gpu_bab_halfspace_verify(net, lb, ub, struct('Hg', []), opts);                      % empty Hg
assert(strcmp(vb, 'skip'), 'empty Hg must skip, got %s', vb);
[vc, ~] = gpu_bab_halfspace_verify(net, lb, ub, struct('Hg', HalfSpace([1 0 0], 1)), opts);   % wrong dim (3 != nOut 2)
assert(strcmp(vc, 'skip'), 'dimension-mismatched Hg must skip, got %s', vc);

%% Test 4: unsupported architecture (tanh) -> 'skip'
ulayers = [featureInputLayer(4,'Name','in'); fullyConnectedLayer(5,'Name','f1'); ...
           tanhLayer('Name','t'); fullyConnectedLayer(2,'Name','out')];
unet = matlab2nnv(dlnetwork(ulayers));
[v4, ~] = gpu_bab_halfspace_verify(unet, lb, ub, struct('Hg', HalfSpace([1 0], 0)), opts);
assert(strcmp(v4, 'skip'), 'unsupported arch (tanh) must skip, got %s', v4);

%% Test 5: a multi-disjunct UNION the reachable set enters -> must NOT be 'robust'
% union = {y_1 <= med1} OR {y_1 >= max1 + 100}: the first is entered, so the union is unsafe.
Hu = [HalfSpace([1 0], med1), HalfSpace([-1 0], -(max1 + 100))];
[v5, ~] = gpu_bab_halfspace_verify(net, lb, ub, struct('Hg', Hu), opts);
assert(~strcmp(v5, 'robust'), 'CARDINAL: robust on a union that the reachable set enters (-150)');

disp('test_soundness_gpu_bab_halfspace: all sections passed');
