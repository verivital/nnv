%% Run examples to load a CNN from different formats 
% As of October 2nd, 2019, only cnn1 is complete. Working on adding support
% to the other layers where the remaining test cases fail.
% Test 1 - matlab
cnn1 = load_cnn('matlab','models/TEST_NET.mat'); % Works
% Test 2 - keras
cnn2 = load_cnn('keras','models/final_model.h5'); % Fails because it has
% a type of layer that is not supported by nnv yet
% Test 3 - onnx
cnn3 = load_cnn('onnx','models/cifarResNet.onnx','classification'); % Fails 
% because it has a type of layer that is not supported by nnv yet
% Test 4 - caffe
cnn4 = load_cnn('caffe','models/deploy_alexnet.prototxt','models/bvlc_alexnet.caffemodel');
% Test 5 - onnx
classes = ["airplane" "automobile" "bird" "cat" "deer" "dog" "frog" "horse" "ship" "truck"];
cnn5 = load_cnn('onnx','models/cifarResNet.onnx','classification',classes);