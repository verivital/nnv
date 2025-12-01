% test_regression_mnist_fc
% Regression test for MNIST fully-connected network verification
% Tests L-infinity adversarial robustness verification
% Based on: examples/Tutorial/NN/MNIST/verify_fc.m
% To run: results = runtests('test_regression_mnist_fc')

%% Test 1: Load MNIST FC network
rng(42);

% Get path to MNIST model
mnist_path = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'Tutorial', filesep, ...
    'NN', filesep, 'MNIST', filesep, 'mnist_model_fc.mat'];

% Load network
mnist_model = load(mnist_path);

assert(isfield(mnist_model, 'net'), 'Model should contain net field');

%% Test 2: Convert to NNV format
rng(42);

mnist_path = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'Tutorial', filesep, ...
    'NN', filesep, 'MNIST', filesep, 'mnist_model_fc.mat'];
mnist_model = load(mnist_path);

% Convert to NNV
net = matlab2nnv(mnist_model.net);

assert(~isempty(net), 'Network should be created');

%% Test 3: Network evaluation
rng(42);

mnist_path = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'Tutorial', filesep, ...
    'NN', filesep, 'MNIST', filesep, 'mnist_model_fc.mat'];
mnist_model = load(mnist_path);
net = matlab2nnv(mnist_model.net);

% Create a simple test image (28x28 grayscale)
test_img = single(128 * ones(28, 28));

% Evaluate
output = net.evaluate(test_img);

assert(length(output) == 10, 'Output should have 10 classes');

%% Test 4: ImageStar input set creation
rng(42);

% Create a test image
img = single(128 * ones(28, 28));

% Create input set with L_inf perturbation
disturbance = 1;  % 1 pixel value disturbance
ones_ = ones(size(img), 'single');

% Create ImageStar with lower and upper bounds (clipped to valid range)
lb_clip = max(img - disturbance, 0);
ub_clip = min(img + disturbance, 255);
IS = ImageStar(lb_clip, ub_clip);

assert(~isempty(IS), 'ImageStar should be created');

%% Test 5: Robustness verification with approx-star
rng(42);

mnist_path = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'Tutorial', filesep, ...
    'NN', filesep, 'MNIST', filesep, 'mnist_model_fc.mat'];
mnist_model = load(mnist_path);
net = matlab2nnv(mnist_model.net);

% Create a test image (all 128s for simplicity)
img = single(128 * ones(28, 28));

% Evaluate to get predicted class
output = net.evaluate(img);
[~, predicted_class] = max(output);
target = single(predicted_class);

% Create input set with very small perturbation for quick verification
disturbance = 0.5;
lb_clip = max(img - disturbance, 0);
ub_clip = min(img + disturbance, 255);
IS = ImageStar(lb_clip, ub_clip);

% Define reachability options
reachOptions = struct;
reachOptions.reachMethod = 'approx-star';

% Verify robustness
res = net.verify_robustness(IS, reachOptions, target);

% Result should be 0 (unknown), 1 (robust), or 2 (not robust)
assert(res >= 0 && res <= 2, 'Verification result should be valid');

%% Test 6: Reachability analysis with approx-star
rng(42);

mnist_path = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'Tutorial', filesep, ...
    'NN', filesep, 'MNIST', filesep, 'mnist_model_fc.mat'];
mnist_model = load(mnist_path);
net = matlab2nnv(mnist_model.net);

% Create a simple input set
img = single(128 * ones(28, 28));
disturbance = 0.5;
lb_clip = max(img - disturbance, 0);
ub_clip = min(img + disturbance, 255);
IS = ImageStar(lb_clip, ub_clip);

% Reachability options
reachOptions = struct;
reachOptions.reachMethod = 'approx-star';

% Compute reachable set
R = net.reach(IS, reachOptions);

assert(~isempty(R), 'Reachable set should be computed');

%% Test 7: Output range computation
rng(42);

mnist_path = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'Tutorial', filesep, ...
    'NN', filesep, 'MNIST', filesep, 'mnist_model_fc.mat'];
mnist_model = load(mnist_path);
net = matlab2nnv(mnist_model.net);

% Create input set
img = single(128 * ones(28, 28));
disturbance = 0.5;
lb_clip = max(img - disturbance, 0);
ub_clip = min(img + disturbance, 255);
IS = ImageStar(lb_clip, ub_clip);

reachOptions = struct;
reachOptions.reachMethod = 'approx-star';
R = net.reach(IS, reachOptions);

% Get output ranges
R_out = net.reachSet{end};
[lb_out, ub_out] = R_out.getRanges;
lb_out = squeeze(lb_out);
ub_out = squeeze(ub_out);

assert(length(lb_out) == 10, 'Should have 10 lower bounds');
assert(length(ub_out) == 10, 'Should have 10 upper bounds');
assert(all(lb_out <= ub_out), 'Lower bounds should be <= upper bounds');

%% Test 8: Soundness check - concrete output within reachable bounds
rng(42);

mnist_path = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'Tutorial', filesep, ...
    'NN', filesep, 'MNIST', filesep, 'mnist_model_fc.mat'];
mnist_model = load(mnist_path);
net = matlab2nnv(mnist_model.net);

% Create input set
img = single(128 * ones(28, 28));
disturbance = 1;
lb_clip = max(img - disturbance, 0);
ub_clip = min(img + disturbance, 255);
IS = ImageStar(lb_clip, ub_clip);

reachOptions = struct;
reachOptions.reachMethod = 'approx-star';
R = net.reach(IS, reachOptions);

% Get output ranges
R_out = net.reachSet{end};
[lb_out, ub_out] = R_out.getRanges;
lb_out = squeeze(lb_out);
ub_out = squeeze(ub_out);

% Evaluate center point and perturbed points
center_output = net.evaluate(img);
corner1_output = net.evaluate(lb_clip);
corner2_output = net.evaluate(ub_clip);

% Check soundness: all concrete outputs should be within bounds
tol = 1e-5;
assert(all(center_output >= lb_out - tol) && all(center_output <= ub_out + tol), ...
    'Center output should be within reachable bounds');
assert(all(corner1_output >= lb_out - tol) && all(corner1_output <= ub_out + tol), ...
    'Corner 1 output should be within reachable bounds');
assert(all(corner2_output >= lb_out - tol) && all(corner2_output <= ub_out + tol), ...
    'Corner 2 output should be within reachable bounds');

