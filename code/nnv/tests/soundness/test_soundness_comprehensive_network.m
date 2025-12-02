% test_soundness_comprehensive_network
% Comprehensive soundness test using multiple layer types
% Tests layer-by-layer soundness verification on a multi-layer network
% To run: results = runtests('test_soundness_comprehensive_network')

%% Test 1: Linear layers only (Conv -> AvgPool -> Flatten -> FC)
% This should pass since all layers are linear
rng(42);

fprintf('\n=== Test 1: Linear Layers Network ===\n');

% Create layers
W_conv = randn(3, 3, 1, 2) * 0.3;  % 3x3 conv, 1->2 channels
b_conv = reshape([0.1, -0.1], [1, 1, 2]);  % Bias must be [1, 1, NumFilters]
L1 = Conv2DLayer(W_conv, b_conv);

L2 = AveragePooling2DLayer([2 2], [2 2], [0 0 0 0]);

L3 = FlattenLayer('flatten');
L3.Type = 'nnet.cnn.layer.FlattenLayer';

W_fc = randn(2, 8) * 0.3;  % 8 inputs (2x2x2), 2 outputs
b_fc = [0.1; -0.1];
L4 = FullyConnectedLayer(W_fc, b_fc);

% Build network
net1 = NN({L1, L2, L3, L4});

% Create input ImageStar (6x6x1 with 2 predicate variables)
V1 = zeros(6, 6, 1, 3);
V1(:,:,1,1) = 0.5 + 0.1*randn(6, 6);  % center
V1(2,2,1,2) = 0.1;  % basis 1
V1(4,4,1,3) = 0.1;  % basis 2

C1 = [eye(2); -eye(2)];
d1 = ones(4, 1);
pred_lb1 = [-1; -1];
pred_ub1 = [1; 1];

input_is1 = ImageStar(V1, C1, d1, pred_lb1, pred_ub1);

% Run layer-by-layer verification
results1 = soundness_test_utils.verify_network_soundness_layer_by_layer(net1, input_is1, 'approx-star', 30, 1e-5, true);

% Assert all layers pass
for i = 1:length(results1)
    assert(results1(i).passed, 'Linear network layer %d (%s) failed soundness', i, results1(i).layer_name);
end
fprintf('Test 1 PASSED: All linear layers are sound.\n');

%% Test 2: Network with ReLU activation
rng(42);

fprintf('\n=== Test 2: Network with ReLU ===\n');

W_conv2 = randn(2, 2, 1, 1) * 0.3;
b_conv2 = reshape(0.1, [1, 1, 1]);
L2_1 = Conv2DLayer(W_conv2, b_conv2);

L2_2 = ReluLayer();

L2_3 = FlattenLayer('flatten2');
L2_3.Type = 'nnet.cnn.layer.FlattenLayer';

W_fc2 = randn(2, 9) * 0.3;  % 3x3 from 4x4->(2x2 conv)->3x3
b_fc2 = [0.1; -0.1];
L2_4 = FullyConnectedLayer(W_fc2, b_fc2);

net2 = NN({L2_1, L2_2, L2_3, L2_4});

% 4x4 input
V2 = zeros(4, 4, 1, 2);
V2(:,:,1,1) = rand(4, 4);
V2(2,2,1,2) = 0.2;

input_is2 = ImageStar(V2, [1; -1], [0.5; 0.5], -0.5, 0.5);

results2 = soundness_test_utils.verify_network_soundness_layer_by_layer(net2, input_is2, 'approx-star', 30, 1e-5, true);

for i = 1:length(results2)
    assert(results2(i).passed, 'ReLU network layer %d (%s) failed soundness', i, results2(i).layer_name);
end
fprintf('Test 2 PASSED: ReLU network is sound.\n');

%% Test 3: Network with Sigmoid activation
rng(42);

fprintf('\n=== Test 3: Network with Sigmoid ===\n');

W_conv3 = randn(2, 2, 1, 1) * 0.3;
b_conv3 = reshape(0, [1, 1, 1]);
L3_1 = Conv2DLayer(W_conv3, b_conv3);

L3_2 = SigmoidLayer();

L3_3 = FlattenLayer('flatten3');
L3_3.Type = 'nnet.cnn.layer.FlattenLayer';

W_fc3 = randn(2, 9) * 0.3;
b_fc3 = [0; 0];
L3_4 = FullyConnectedLayer(W_fc3, b_fc3);

net3 = NN({L3_1, L3_2, L3_3, L3_4});

V3 = zeros(4, 4, 1, 2);
V3(:,:,1,1) = 0.5*ones(4, 4);  % Center near sigmoid's linear region
V3(2,2,1,2) = 0.1;

input_is3 = ImageStar(V3, [1; -1], [0.3; 0.3], -0.3, 0.3);

results3 = soundness_test_utils.verify_network_soundness_layer_by_layer(net3, input_is3, 'approx-star', 30, 1e-5, true);

for i = 1:length(results3)
    assert(results3(i).passed, 'Sigmoid network layer %d (%s) failed soundness', i, results3(i).layer_name);
end
fprintf('Test 3 PASSED: Sigmoid network is sound.\n');

%% Test 4: Network with Tanh activation
rng(42);

fprintf('\n=== Test 4: Network with Tanh ===\n');

W_conv4 = randn(2, 2, 1, 1) * 0.3;
b_conv4 = reshape(0, [1, 1, 1]);
L4_1 = Conv2DLayer(W_conv4, b_conv4);

L4_2 = TanhLayer();

L4_3 = FlattenLayer('flatten4');
L4_3.Type = 'nnet.cnn.layer.FlattenLayer';

W_fc4 = randn(2, 9) * 0.3;
b_fc4 = [0; 0];
L4_4 = FullyConnectedLayer(W_fc4, b_fc4);

net4 = NN({L4_1, L4_2, L4_3, L4_4});

V4 = zeros(4, 4, 1, 2);
V4(:,:,1,1) = zeros(4, 4);  % Center at 0 for tanh
V4(2,2,1,2) = 0.1;

input_is4 = ImageStar(V4, [1; -1], [0.3; 0.3], -0.3, 0.3);

results4 = soundness_test_utils.verify_network_soundness_layer_by_layer(net4, input_is4, 'approx-star', 30, 1e-5, true);

for i = 1:length(results4)
    assert(results4(i).passed, 'Tanh network layer %d (%s) failed soundness', i, results4(i).layer_name);
end
fprintf('Test 4 PASSED: Tanh network is sound.\n');

%% Test 5: Network with LeakyReLU
rng(42);

fprintf('\n=== Test 5: Network with LeakyReLU ===\n');

W_conv5 = randn(2, 2, 1, 1) * 0.3;
b_conv5 = reshape(0, [1, 1, 1]);
L5_1 = Conv2DLayer(W_conv5, b_conv5);

L5_2 = LeakyReluLayer('leakyrelu', 1, {'in'}, 1, {'out'}, 0.01);

L5_3 = FlattenLayer('flatten5');
L5_3.Type = 'nnet.cnn.layer.FlattenLayer';

W_fc5 = randn(2, 9) * 0.3;
b_fc5 = [0; 0];
L5_4 = FullyConnectedLayer(W_fc5, b_fc5);

net5 = NN({L5_1, L5_2, L5_3, L5_4});

V5 = zeros(4, 4, 1, 2);
V5(:,:,1,1) = rand(4, 4) - 0.5;
V5(2,2,1,2) = 0.2;

input_is5 = ImageStar(V5, [1; -1], [0.5; 0.5], -0.5, 0.5);

results5 = soundness_test_utils.verify_network_soundness_layer_by_layer(net5, input_is5, 'approx-star', 30, 1e-5, true);

for i = 1:length(results5)
    assert(results5(i).passed, 'LeakyReLU network layer %d (%s) failed soundness', i, results5(i).layer_name);
end
fprintf('Test 5 PASSED: LeakyReLU network is sound.\n');

%% Test 6: Network with MaxPooling (exact-star)
rng(42);

fprintf('\n=== Test 6: Network with MaxPooling ===\n');

W_conv6 = randn(2, 2, 1, 1) * 0.3;
b_conv6 = reshape(0.1, [1, 1, 1]);
L6_1 = Conv2DLayer(W_conv6, b_conv6);

L6_2 = MaxPooling2DLayer([2 2], [2 2], [0 0 0 0]);

L6_3 = FlattenLayer('flatten6');
L6_3.Type = 'nnet.cnn.layer.FlattenLayer';

W_fc6 = randn(2, 4) * 0.3;  % 2x2 after pooling
b_fc6 = [0; 0];
L6_4 = FullyConnectedLayer(W_fc6, b_fc6);

net6 = NN({L6_1, L6_2, L6_3, L6_4});

% Small input for MaxPool
V6 = zeros(5, 5, 1, 2);
V6(:,:,1,1) = rand(5, 5);
V6(2,2,1,2) = 0.1;

input_is6 = ImageStar(V6, [1; -1], [0.3; 0.3], -0.3, 0.3);

% MaxPooling needs exact-star for soundness
results6 = soundness_test_utils.verify_network_soundness_layer_by_layer(net6, input_is6, 'exact-star', 30, 1e-5, true);

for i = 1:length(results6)
    assert(results6(i).passed, 'MaxPool network layer %d (%s) failed soundness', i, results6(i).layer_name);
end
fprintf('Test 6 PASSED: MaxPooling network is sound.\n');

%% Test 7: Deep network (Conv -> ReLU -> Conv -> ReLU -> Pool -> Flatten -> FC)
rng(42);

fprintf('\n=== Test 7: Deep Network ===\n');

W_c1 = randn(3, 3, 1, 2) * 0.2;
b_c1 = reshape([0.1, -0.1], [1, 1, 2]);
L7_1 = Conv2DLayer(W_c1, b_c1);

L7_2 = ReluLayer();

W_c2 = randn(2, 2, 2, 1) * 0.2;
b_c2 = reshape(0, [1, 1, 1]);
L7_3 = Conv2DLayer(W_c2, b_c2);

L7_4 = ReluLayer();

L7_5 = AveragePooling2DLayer([2 2], [2 2], [0 0 0 0]);

L7_6 = FlattenLayer('flatten7');
L7_6.Type = 'nnet.cnn.layer.FlattenLayer';

W_fc7 = randn(3, 4) * 0.2;  % Depends on input size
b_fc7 = zeros(3, 1);
L7_7 = FullyConnectedLayer(W_fc7, b_fc7);

net7 = NN({L7_1, L7_2, L7_3, L7_4, L7_5, L7_6, L7_7});

% 8x8 input to get proper output sizes
V7 = zeros(8, 8, 1, 2);
V7(:,:,1,1) = 0.5 + 0.1*randn(8, 8);
V7(3,3,1,2) = 0.1;

input_is7 = ImageStar(V7, [1; -1], [0.5; 0.5], -0.5, 0.5);

results7 = soundness_test_utils.verify_network_soundness_layer_by_layer(net7, input_is7, 'approx-star', 30, 1e-5, true);

for i = 1:length(results7)
    assert(results7(i).passed, 'Deep network layer %d (%s) failed soundness', i, results7(i).layer_name);
end
fprintf('Test 7 PASSED: Deep network is sound.\n');

%% Summary
fprintf('\n=== All Comprehensive Network Soundness Tests PASSED ===\n');
