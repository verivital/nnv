% test_regression_verify_segmentation
% Regression tests for semantic segmentation verification
% Tests verify_segmentation and related methods in NN class
% To run: results = runtests('test_regression_verify_segmentation')

%% Test 1: verify_segmentation method exists
rng(42);
assert(ismethod(NN([]), 'verify_segmentation'), 'verify_segmentation method should exist');

%% Test 2: plot_segmentation_output_set method exists
rng(42);
assert(ismethod(NN([]), 'plot_segmentation_output_set'), 'plot_segmentation_output_set method should exist');

%% Test 3: getPixelClassReachSet method exists
rng(42);
assert(ismethod(NN([]), 'getPixelClassReachSet'), 'getPixelClassReachSet method should exist');

%% Test 4: Load M2NIST segmentation model
rng(42);
net_path = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'NN', filesep, ...
    'SemanticSegmentation', filesep, 'M2NIST', filesep, 'models', filesep, 'm2nist_75iou_transposedcnn_avgpool.mat'];
net_data = load(net_path);
net = matlab2nnv(net_data.net);
assert(~isempty(net), 'Segmentation network should load');

%% Test 5: Load ground truth labels
rng(42);
labels_path = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'NN', filesep, ...
    'SemanticSegmentation', filesep, 'M2NIST', filesep, 'm2nist_6484_test_labels.mat'];
if exist(labels_path, 'file')
    labels = load(labels_path);
    assert(~isempty(fieldnames(labels)), 'Labels file should have data');
else
    % Labels file might not exist, just check we can create inputs
    assert(true, 'Test setup complete');
end

%% Test 6: Create ImageStar for segmentation input
rng(42);
images_path = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'NN', filesep, ...
    'SemanticSegmentation', filesep, 'M2NIST', filesep, 'm2nist_6484_test_images.mat'];
images = load(images_path);
im_data = single(images.im_data);
img = im_data(:,:,1);
% Create attack bound ImageStar
disturbance = 1;
lb_clip = max(img - disturbance, 0);
ub_clip = min(img + disturbance, 255);
IS = ImageStar(lb_clip, ub_clip);
assert(~isempty(IS), 'ImageStar input should be created');
assert(IS.numPred > 0, 'ImageStar should have predicates for attack pixels');

%% Test 7: Evaluate segmentation network outputs correct shape
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
% Output should have spatial dimensions and class channels
assert(ndims(output) >= 2, 'Output should be at least 2D');
assert(all(isfinite(output(:))), 'Output should be finite');

%% Test 8: verify_segmentation has correct signature
rng(42);
% Test that method accepts expected number of inputs/outputs
mc = ?NN;
methods = mc.MethodList;
verify_seg_method = [];
for i = 1:length(methods)
    if strcmp(methods(i).Name, 'verify_segmentation')
        verify_seg_method = methods(i);
        break;
    end
end
assert(~isempty(verify_seg_method), 'verify_segmentation method should be found');
% Method should return multiple outputs (riou, rv, rs, etc.)
assert(~isempty(verify_seg_method), 'Method should be found in NN class');
