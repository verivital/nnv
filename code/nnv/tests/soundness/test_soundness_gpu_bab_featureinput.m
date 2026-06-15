% test_soundness_gpu_bab_featureinput
% nn_to_ops must accept a leading FeatureInputLayer (flat [n] feature input, e.g. the FC nets
% in sat_relu / safenlp imported with InputDataFormats="BC") as an engine-input passthrough,
% and the resulting op-list must compute the SAME function as net.evaluate -- the orientation
% guard's soundness invariant (gpu_bab_ibp on a degenerate box == net.evaluate). This is the
% import-correctness check that gates the GPU-BaB FC unsat pre-check on these benchmarks.

%% Test 1: FeatureInputLayer FC+ReLU net -> nn_to_ops op-list == net.evaluate (degenerate IBP)
rng(1);
layers = [featureInputLayer(8,'Name','in')
          fullyConnectedLayer(16,'Name','fc1'); reluLayer('Name','r1')
          fullyConnectedLayer(12,'Name','fc2'); reluLayer('Name','r2')
          fullyConnectedLayer(5,'Name','fc3')];
net = dlnetwork(layers);
nnvnet = matlab2nnv(net);
assert(any(cellfun(@(L) contains(class(L),'FeatureInput'), nnvnet.Layers)), ...
    'expected a FeatureInputLayer in the nnv net');
ops = nn_to_ops(nnvnet, 'colmajor');         % must not error on FeatureInputLayer
maxe = 0;
for s = 1:8
    x = randn(8,1);
    yo = gpu_bab_ibp(ops, x, x, 'double'); yo = yo(:);
    yn = nnvnet.evaluate(x); yn = yn(:);
    maxe = max(maxe, max(abs(yo - yn)));
end
assert(maxe < 1e-6, sprintf('op-list (FeatureInput) must match net.evaluate; max err %.3e', maxe));

%% Test 2: a flat normaffine (the ElementwiseAffine / flat-BN fold) is exact in the op-list
% Build an op-list by hand: affine -> flat normaffine -> relu -> affine, and check the
% normaffine (y = s.*x + t) is applied exactly by the IBP/degenerate eval.
rng(2);
W1 = randn(6,4); b1 = randn(6,1); s = randn(6,1); t = randn(6,1); W2 = randn(3,6); b2 = randn(3,1);
ops = {struct('type','affine','W',W1,'b',b1,'src',0), ...
       struct('type','normaffine','scale',s,'shift',t,'shape',[6 1 1],'src',1), ...
       struct('type','relu','src',2), ...
       struct('type','affine','W',W2,'b',b2,'src',3)};
maxe = 0;
for k = 1:8
    x = randn(4,1);
    ref = W2 * max(s.*(W1*x+b1)+t, 0) + b2;          % the exact composed function
    yo = gpu_bab_ibp(ops, x, x, 'double'); yo = yo(:);
    maxe = max(maxe, max(abs(yo - ref)));
end
assert(maxe < 1e-9, sprintf('flat normaffine must be exact in the op-list; max err %.3e', maxe));

%% Summary
disp('test_soundness_gpu_bab_featureinput: all sections passed');
