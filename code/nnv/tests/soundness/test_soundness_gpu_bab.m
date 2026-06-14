% test_soundness_gpu_bab
% Soundness + correctness tests for the GPU-BaB engine (engine/nn/gpu_bab).
% Verifies that the LP-free bound-propagation + branch-and-bound cores are SOUND
% (over-approximate / sound-or-unknown) and that the tightenings actually tighten.
% To run: results = runtests('test_soundness_gpu_bab')
%
% NOTE: each %% section runs independently under runtests, so every quantity shared
% across sections is computed here in the shared-variables setup (not in a prior test).

% ---- shared setup (deterministic): fixed FC+ReLU net, box, MC outputs, baseline bounds ----
% (engine/nn/gpu_bab is on the path via startup_nnv / install)
rng(0);
W1 = randn(16,5);  b1 = randn(16,1);
W2 = randn(16,16); b2 = randn(16,1);
W3 = randn(3,16);  b3 = randn(3,1);
net = NN({FullyConnectedLayer(W1,b1), ReluLayer(), ...
          FullyConnectedLayer(W2,b2), ReluLayer(), ...
          FullyConnectedLayer(W3,b3)});
ops = nn_to_ops(net);
lb  = -ones(5,1); ub = ones(5,1);
Cspec = [1 -1 0; 1 0 -1];                       % robustness spec for class 1 (e_1 - e_j)
[ilb, iub] = gpu_bab_ibp(ops, lb, ub, 'double');% IBP output box
Xs = lb + (ub - lb) .* rand(5, 5000);           % Monte-Carlo samples
Ys = gpu_bab_ibp(ops, Xs, Xs, 'double');        % exact outputs (degenerate box = exact)
trueMin = min(Cspec * Ys, [], 2);               % empirical min of the spec over the box
m1 = gpu_bab_crown_spec(ops, lb, ub, Cspec, 'double');  % fixed-slope CROWN spec margin

%% Test 1: op extraction reproduces NNV evaluate exactly
x = lb + (ub - lb) .* rand(5,1);
yo = gpu_bab_ibp(ops, x, x, 'double');
yn = net.evaluate(x);
assert(max(abs(yo - yn(:))) < 1e-9, 'nn_to_ops + gpu_bab_ibp must match NNV evaluate');

%% Test 2: IBP is sound (output box contains every sampled output)
assert(all(Ys >= ilb - 1e-6 & Ys <= iub + 1e-6, 'all'), 'IBP box must contain all sampled outputs');

%% Test 3: CROWN is sound and no looser than IBP
[clb, cub] = gpu_bab_crown(ops, lb, ub, 'double');
assert(all(Ys >= clb - 1e-6 & Ys <= cub + 1e-6, 'all'), 'CROWN box must contain all sampled outputs');
assert(mean(cub - clb) <= mean(iub - ilb) + 1e-9, 'CROWN must be no looser than IBP');

%% Test 4: batched CROWN spec == single-box, and is a sound lower bound
m4 = gpu_bab_crown_spec(ops, [lb lb lb lb], [ub ub ub ub], Cspec, 'double');
assert(max(abs(m4 - m1), [], 'all') < 1e-9, 'batched spec must equal single-box per column');
assert(all(m1 <= trueMin + 1e-4), 'spec margin must lower-bound the true min over the box');

%% Test 5: CROWN-tight is sound and no looser than IBP-intermediate CROWN
mt = gpu_bab_crown_tight(ops, lb, ub, Cspec, 'double');
assert(all(mt <= trueMin + 1e-4), 'CROWN-tight margin must be sound (<= true min)');
assert(min(mt) >= min(m1) - 1e-6, 'CROWN-tight must be no looser than IBP-intermediate CROWN');

%% Test 6: alpha-CROWN tightens (>= fixed-slope) and stays sound
ma = gpu_bab_crown_alpha(ops, lb, ub, Cspec, 'double', 30, 0.3);
assert(all(ma >= m1 - 1e-4), 'alpha-CROWN must be no looser than fixed-slope CROWN');
assert(all(ma <= trueMin + 1e-4), 'alpha-CROWN margin must be sound (<= true min)');

%% Test 7: ReLU-split robust verdict is sound (no misclassifying sample)
x0 = lb + (ub - lb) .* rand(5,1);
[~, tl] = max(net.evaluate(x0));
opts = struct('precision','double', 'maxNodes',500, 'margin',1e-4);
s = gpu_bab_relu_split(ops, x0 - 0.01, x0 + 0.01, tl, 3, opts);
assert(ismember(s, {'robust','unsafe','unknown'}), 'verdict must be robust/unsafe/unknown');
Xr = (x0 - 0.01) + 0.02 .* rand(5, 8000);
Yr = gpu_bab_ibp(ops, Xr, Xr, 'double');
[~, pr] = max(Yr, [], 1);
assert(~strcmp(s, 'robust') || all(pr == tl), 'a robust verdict must have no misclassifying sample');

%% Test 8: ReLU-split unsafe verdict carries a real counterexample
x0b = lb + (ub - lb) .* rand(5,1);
[~, tlb] = max(net.evaluate(x0b));
opts2 = struct('precision','double', 'maxNodes',500, 'margin',1e-4, 'nSample',60);
[s2, info2] = gpu_bab_relu_split(ops, x0b - 5, x0b + 5, tlb, 3, opts2);
ok8 = ~strcmp(s2, 'unsafe');
if strcmp(s2, 'unsafe')
    yc = net.evaluate(info2.cex);
    [~, pc] = max(yc);
    ok8 = (pc ~= tlb);
end
assert(ok8, 'an unsafe verdict must carry an input that actually misclassifies');

%% Test 9: JOINT alpha-CROWN is sound and no looser than fixed-slope CROWN-tight
% gpu_bab_crown_alpha_joint optimizes the unstable-ReLU lower slopes INSIDE every
% intermediate backward pass (not just the final spec), initialized at the min-area
% slopes so iter 0 reproduces gpu_bab_crown_tight, then ascends + keeps-best.
mt9 = gpu_bab_crown_tight(ops, lb, ub, Cspec, 'double');
[mj9, ij9] = gpu_bab_crown_alpha_joint(ops, lb, ub, Cspec, 'double', 30, 0.5);
assert(all(mj9 <= trueMin + 1e-4), 'joint alpha-CROWN margin must be sound (<= true min)');
assert(ij9.alpha_minmargin >= sum(min(mt9,[],1)) - 1e-6, 'joint alpha-CROWN must be no looser than fixed-slope CROWN-tight');

%% Summary
disp('test_soundness_gpu_bab: all sections passed');
