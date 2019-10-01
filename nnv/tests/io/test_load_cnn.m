%% Run examples to load a CNN from different formats 
% Test 1 - matlab
cnn1 = load_cnn('matlab','models/TEST_NET.mat');
% Test 2 - keras
cnn2 = load_cnn('keras','models/final_model.h5');
% Test 3 - onnx
cnn3 = load_cnn('onnx','models/cifarResNet.onnx','classification');
% Test 4 - caffe
cnn4 = load_cnn('caffe','models/deploy_alexnet.prototxt','models/bvlc_alexnet.caffemodel');
% Test 5 - onnx
classes = ["airplane" "automobile" "bird" "cat" "deer" "dog" "frog" "horse" "ship" "truck"];
cnn5 = load_cnn('onnx','models/cifarResNet.onnx','classification',classes);