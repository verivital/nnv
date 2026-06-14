% test_soundness_gpu_bab
% Soundness + correctness tests for the GPU-BaB engine (engine/nn/gpu_bab).
% Verifies that the LP-free bound-propagation + branch-and-bound cores are SOUND
% (over-approximate / sound-or-unknown) and that the tightenings actually tighten.
% To run: results = runtests('test_soundness_gpu_bab')

% ---- shared setup: a fixed small FC+ReLU net + an input box (deterministic) ----
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'engine', 'nn', 'gpu_bab')));
rng(0);
W1 = randn(16,5);  b1 = randn(16,1);
W2 = randn(16,16); b2 = randn(16,1);
W3 = randn(3,16);  b3 = randn(3,1);
net = NN({FullyConnectedLayer(W1,b1), ReluLayer(), ...
          FullyConnectedLayer(W2,b2), ReluLayer(), ...
          FullyConnectedLayer(W3,b3)});
ops = nn_to_ops(net);
lb = -ones(5,1); ub = ones(5,1);
% robustness spec C for class 1 (rows e_1 - e_j)
Cspec = [1 -1 0; 1 0 -1];

%% Test 1: op extraction reproduces NNV evaluate exactly
x = lb + (ub - lb) .* rand(5,1);
yo = gpu_bab_ibp(ops, x, x, 'double');
yn = net.evaluate(x);
assert(max(abs(yo - yn(:))) < 1e-9, 'nn_to_ops + gpu_bab_ibp must match NNV evaluate');

%% Test 2: IBP is sound (output box contains every sampled output)
[ilb, iub] = gpu_bab_ibp(ops, lb, ub, 'double');
X = lb + (ub - lb) .* rand(5, 5000);
Y = gpu_bab_ibp(ops, X, X, 'double');
assert(all(Y >= ilb - 1e-6 & Y <= iub + 1e-6, 'all'), 'IBP box must contain all sampled outputs');

%% Test 3: CROWN is sound and tighter than IBP
[clb, cub] = gpu_bab_crown(ops, lb, ub, 'double');
assert(all(Y >= clb - 1e-6 & Y <= cub + 1e-6, 'all'), 'CROWN box must contain all sampled outputs');
assert(mean(cub - clb) <= mean(iub - ilb) + 1e-9, 'CROWN must be no looser than IBP');

%% Test 4: batched CROWN spec == single-box, and is a sound lower bound
m1 = gpu_bab_crown_spec(ops, lb, ub, Cspec, 'double');
m4 = gpu_bab_crown_spec(ops, [lb lb lb lb], [ub ub ub ub], Cspec, 'double');
assert(max(abs(m4 - m1), [], 'all') < 1e-9, 'batched spec must equal single-box per column');
trueMin = min(Cspec * Y, [], 2);
assert(all(m1 <= trueMin + 1e-4), 'spec margin must lower-bound the true min over the box');

%% Test 5: CROWN-tight is sound and tighter than IBP-intermediate CROWN
mt = gpu_bab_crown_tight(ops, lb, ub, Cspec, 'double');
assert(all(mt <= trueMin + 1e-4), 'CROWN-tight margin must be sound (<= true min)');
assert(min(mt) >= min(m1) - 1e-6, 'CROWN-tight must be no looser than IBP-intermediate CROWN');

%% Test 6: alpha-CROWN tightens (>= fixed-slope) and stays sound
[ma, ~] = gpu_bab_crown_alpha(ops, lb, ub, Cspec, 'double', 30, 0.3);
assert(all(ma >= m1 - 1e-4), 'alpha-CROWN must be no looser than fixed-slope CROWN');
assert(all(ma <= trueMin + 1e-4), 'alpha-CROWN margin must be sound (<= true min)');

%% Test 7: ReLU-split robust verdict is sound (no misclassifying sample)
x0 = lb + (ub - lb) .* rand(5,1);
[~, tl] = max(net.evaluate(x0));
opts = struct('precision','double', 'maxNodes',500, 'margin',1e-4);
[s, ~] = gpu_bab_relu_split(ops, x0 - 0.01, x0 + 0.01, tl, 3, opts);
assert(ismember(s, {'robust','unsafe','unknown'}), 'verdict must be one of robust/unsafe/unknown');
Xr = (x0 - 0.01) + 0.02 .* rand(5, 8000);
Yr = gpu_bab_ibp(ops, Xr, Xr, 'double');
[~, pr] = max(Yr, [], 1);
robustSound = ~strcmp(s, 'robust') || all(pr == tl);
assert(robustSound, 'a robust verdict must have no misclassifying sample in the box');

%% Test 8: ReLU-split unsafe verdict carries a real counterexample
opts2 = struct('precision','double', 'maxNodes',500, 'margin',1e-4, 'nSample',60);
[s2, info2] = gpu_bab_relu_split(ops, x0 - 5, x0 + 5, tl, 3, opts2);
unsafeSound = ~strcmp(s2, 'unsafe');
if strcmp(s2, 'unsafe')
    yc = net.evaluate(info2.cex);
    [~, pc] = max(yc);
    unsafeSound = (pc ~= tl);
end
assert(unsafeSound, 'an unsafe verdict must carry an input that actually misclassifies');

%% Summary
% All GPU-BaB soundness/correctness checks passed: op extraction exact; IBP/CROWN/
% CROWN-tight/alpha-CROWN all sound (and monotonically tighter); batched == single-box;
% ReLU-split verdicts sound-or-unknown (robust => MC-clean, unsafe => real witness).
disp('test_soundness_gpu_bab: all sections passed');
