% test_soundness_gpu_bab_conv_batched
% gpu_bab_crown_spec_dag (the batched, conv-capable spec bounding) must be SOUND on a
% sequential conv net: the conv/avgpool adjoint is EXACT (a degenerate box -> CROWN == the
% concrete C*evaluate), the bound is a valid lower bound (<= every sampled output), and the
% batched conv ReLU-split never CONTRADICTS the serial gpu_bab_relu_split (tight). Synthetic
% MNIST-style convnet (conv->relu->avgpool->conv->relu->fc). CPU-only so it runs in GPU-less
% CI; loops are in local functions to avoid the %%-section script-test line-number quirk.

%% Test 1: degenerate CROWN-dag == C*evaluate (the batched conv/avgpool adjoint is exact)
[ops, nnvnet] = i_build_convnet();
xt = i_fixed_input();
yn = nnvnet.evaluate(reshape(xt, [8 8 1])); yn = yn(:);
C = [eye(9) -ones(9,1)];
md = gpu_bab_crown_spec_dag(ops, xt(:), xt(:), C, 'double');
assert(max(abs(md(:) - C*yn)) < 1e-5, ...
    'degenerate conv CROWN-dag must equal C*evaluate (the conv/avgpool adjoint must be exact)');

%% Test 2: CROWN-dag is a sound lower bound over a box (<= the Monte-Carlo min)
[ops, nnvnet] = i_build_convnet();
xt = i_fixed_input();
C = [eye(9) -ones(9,1)];
lb = xt(:) - 0.05; ub = xt(:) + 0.05;
trueMin = i_mc_min(nnvnet, C, lb, ub, 1500);
m = gpu_bab_crown_spec_dag(ops, lb, ub, C, 'double');
assert(all(m(:) <= trueMin + 1e-4), 'conv CROWN-dag must be a sound lower bound (<= true min over the box)');

%% Test 3: batched conv ReLU-split never CONTRADICTS serial (sound-or-unknown)
[ops, nnvnet] = i_build_convnet();
xt = i_fixed_input();
yn = nnvnet.evaluate(reshape(xt, [8 8 1])); yn = yn(:); [~, tl] = max(yn);
nContra = i_compare_conv(ops, xt, tl);
assert(nContra == 0, sprintf('batched vs serial conv ReLU-split: %d robust/unsafe contradictions', nContra));

%% Summary
disp('test_soundness_gpu_bab_conv_batched: all sections passed');

% ----------------------------------------------------------------------------------------
function [ops, nnvnet] = i_build_convnet()
% Small sequential MNIST-style convnet (8x8x1 -> conv -> relu -> avgpool -> conv -> relu -> fc10).
    rng(5);
    lg = layerGraph([
        imageInputLayer([8 8 1],'Normalization','none','Name','in')
        convolution2dLayer(3,4,'Padding',1,'Name','c1'); reluLayer('Name','r1')
        averagePooling2dLayer(2,'Stride',2,'Name','p1')
        convolution2dLayer(3,8,'Padding',1,'Name','c2'); reluLayer('Name','r2')
        fullyConnectedLayer(10,'Name','fc')]);
    nnvnet = matlab2nnv(dlnetwork(lg));
    ops = nn_to_ops(nnvnet, 'colmajor');   % MATLAB-native dlnetwork -> column-major flatten
end

function x = i_fixed_input()
    rng(7); x = rand(64, 1);               % 8x8x1 flat input, deterministic
end

function mn = i_mc_min(nnvnet, C, lb, ub, nS)
    mn = inf(size(C,1), 1);
    rng(11); X = lb + (ub - lb) .* rand(numel(lb), nS);
    for s = 1:nS
        ys = nnvnet.evaluate(reshape(X(:,s), [8 8 1]));
        mn = min(mn, C * ys(:));
    end
end

function nContra = i_compare_conv(ops, xt, tl)
% serial(tight) vs batched(DAG) ReLU-split over a few eps; count robust/unsafe contradictions.
    nContra = 0; nC = 10;
    for ep = [0.01 0.03 0.06]
        lb = single(xt(:) - ep); ub = single(xt(:) + ep);
        [sS,~] = gpu_bab_relu_split(ops, lb, ub, tl, nC, struct('precision','single','maxNodes',1500));
        [sB,~] = gpu_bab_relu_split_batched(ops, lb, ub, tl, nC, struct('precision','single','maxNodes',1500,'maxFrontier',128));
        if (strcmp(sS,'robust') && strcmp(sB,'unsafe')) || (strcmp(sS,'unsafe') && strcmp(sB,'robust'))
            nContra = nContra + 1;
        end
    end
end
