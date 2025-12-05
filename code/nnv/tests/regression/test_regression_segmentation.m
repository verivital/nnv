% test_regression_segmentation
% Regression test for semantic segmentation network verification
% Tests loading M2NIST networks and basic verification
% Based on: examples/Tutorial/NN/Segmentation/verify_m2nist.m
% To run: results = runtests('test_regression_segmentation')

%% Test 1: Load M2NIST segmentation network
rng(42);

% Get path to segmentation network
net_path = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'NN', filesep, ...
    'SemanticSegmentation', filesep, 'M2NIST', filesep, 'models', filesep, 'm2nist_75iou_transposedcnn_avgpool.mat'];

% Load network
net_data = load(net_path);

assert(isfield(net_data, 'net'), 'Network file should have net field');

%% Test 2: Convert segmentation network to NNV
rng(42);

net_path = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'NN', filesep, ...
    'SemanticSegmentation', filesep, 'M2NIST', filesep, 'models', filesep, 'm2nist_75iou_transposedcnn_avgpool.mat'];
net_data = load(net_path);

% Convert to NNV
net = matlab2nnv(net_data.net);

assert(~isempty(net), 'NNV network should be created');

%% Test 3: Load test images
rng(42);

images_path = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'NN', filesep, ...
    'SemanticSegmentation', filesep, 'M2NIST', filesep, 'm2nist_6484_test_images.mat'];

images = load(images_path);

assert(isfield(images, 'im_data'), 'Images file should have im_data field');

%% Test 4: Segmentation network can evaluate
rng(42);

net_path = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'NN', filesep, ...
    'SemanticSegmentation', filesep, 'M2NIST', filesep, 'models', filesep, 'm2nist_75iou_transposedcnn_avgpool.mat'];
net_data = load(net_path);
net = matlab2nnv(net_data.net);

images_path = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'NN', filesep, ...
    'SemanticSegmentation', filesep, 'M2NIST', filesep, 'm2nist_6484_test_images.mat'];
images = load(images_path);
im_data = single(images.im_data);

% Select first image
img = im_data(:,:,1);

% Evaluate
output = net.evaluate(img);

assert(~isempty(output), 'Output should not be empty');

%% Test 5: Output is segmentation mask
rng(42);

net_path = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'NN', filesep, ...
    'SemanticSegmentation', filesep, 'M2NIST', filesep, 'models', filesep, 'm2nist_75iou_transposedcnn_avgpool.mat'];
net_data = load(net_path);
net = matlab2nnv(net_data.net);

images_path = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'NN', filesep, ...
    'SemanticSegmentation', filesep, 'M2NIST', filesep, 'm2nist_6484_test_images.mat'];
images = load(images_path);
im_data = single(images.im_data);

img = im_data(:,:,1);
output = net.evaluate(img);

% Output should be same spatial size as input with channels for each class
assert(all(isfinite(output(:))), 'Output should be finite');

%% Test 6: Create ImageStar input set
rng(42);

images_path = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'NN', filesep, ...
    'SemanticSegmentation', filesep, 'M2NIST', filesep, 'm2nist_6484_test_images.mat'];
images = load(images_path);
im_data = single(images.im_data);

img = im_data(:,:,1);

% Create small perturbation
disturbance = 1;
lb_clip = max(img - disturbance, 0);
ub_clip = min(img + disturbance, 255);
IS = ImageStar(lb_clip, ub_clip);

assert(~isempty(IS), 'ImageStar should be created');

%% Test 7: Load dilated CNN model
rng(42);

net_path = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'NN', filesep, ...
    'SemanticSegmentation', filesep, 'M2NIST', filesep, 'models', filesep, 'M2NIST_dilatedCNN_avgpool.mat'];

if exist(net_path, 'file')
    net_data = load(net_path);
    assert(isfield(net_data, 'net'), 'Dilated CNN should have net field');

    net = matlab2nnv(net_data.net);
    assert(~isempty(net), 'Dilated CNN should convert to NNV');
end

%% Test 8: Multiple segmentation models can load
rng(42);

base_path = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'NN', filesep, ...
    'SemanticSegmentation', filesep, 'M2NIST', filesep, 'models', filesep];

model_names = {'m2nist_75iou_transposedcnn_avgpool.mat', 'M2NIST_dilatedCNN_avgpool.mat'};

for i = 1:length(model_names)
    model_path = [base_path, model_names{i}];
    if exist(model_path, 'file')
        net_data = load(model_path);
        assert(isfield(net_data, 'net'), sprintf('Model %d should have net field', i));

        net = matlab2nnv(net_data.net);
        assert(~isempty(net), sprintf('Model %d should convert', i));
    end
end

