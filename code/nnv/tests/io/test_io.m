%to run this as a test, use results_io=runtests('test_io')
%requirements: file must start or end with test
%each test starts with two percent signs followed by the name
%shared vairables must appear before first test
%variables made by a test are not available to other tests.


%shared variables


% Run examples to load a FNN from different formats using nnmt
[v,e,l] = pyversion;
% Check OS running on
% osc = computer;
% if contains(osc,'WIN')
%     sh = '\'; % windows
% else 
%     sh = '/'; % mac and linux
% end
sh = filesep;
% Find the path to the python to use, in case that python does not have all
% dependencies needed, look at help pyversion and py.sys.path to choose the
% correct python environment
e = split(string(e),sh);
e = erase(e,e(end));
pypath = strjoin(e,sh);

% Get paths to the inputs for the first example
genp = split(string(pwd),sh);
rp = {'nnv','code','nnv','tests','io'};

if ~isequal(genp(end-4:end),string(rp)')
  error('Executing from %s. \nPlease, change your current folder path to ../nnv/code/nnv/io/tests and run script again',pwd);
end

fp = genp(1:end-2);
fp = strjoin(fp,filesep);
fp = char(fp);




%___________________________________________________________________________________________________
%tests below originally taken from test_load_cnn.m

%% test 1: load cnn matlab
cnn1 = load_cnn('matlab','models/TEST_NET.mat'); % Works

    
%% test 2: load cnn keras
cnn2 = load_cnn('keras','models/final_model.h5'); % Fails because it
		% has a type of layer that is not supported by nnv yet



%% test 3: load cnn onnx
cnn3 = load_cnn('onnx','models/cifarResNet.onnx','classification'); % Fails because it has a type of layer that is not supported by nnv yet


%% test 4: load cnn caffe
cnn4 = load_cnn('caffe','models/deploy_alexnet.prototxt','models/bvlc_alexnet.caffemodel');

%% test 5: load cnn onnx 2
classes = ["airplane" "automobile" "bird" "cat" "deer" "dog" "frog" "horse" "ship" "truck"];
cnn5 = load_cnn('onnx','models/cifarResNet.onnx','classification', classes);



%___________________________________________________________________________________________________
%tests below originally taken from test_load_nn.m




%% test 6: load nn 1
cont1 = Load_nn([fp sh 'engine' sh 'nnmt' sh 'original_networks' sh 'neural_network_information_13'],[fp sh 'engine' sh 'nnmt' sh 'translated_networks']);

%% test 7: load nn 2
cont2 = Load_nn([fp sh 'engine' sh 'nnmt' sh 'original_networks' sh 'neural_network_information_2'],[fp sh 'engine' sh 'nnmt' sh 'translated_networks']);

%% test 8: load nn 3
cont3 = Load_nn([fp sh 'engine' sh 'nnmt' sh 'original_networks' sh 'MC_16sigX16sigX1tanh.h5'],[fp sh 'engine' sh 'nnmt' sh 'translated_networks']);

%% test 9: load nn 4
cont4 = Load_nn([fp sh 'engine' sh 'nnmt' sh 'original_networks' sh 'ACASXU_run2a_4_3_batch_2000.nnet'],[fp sh 'engine' sh 'nnmt' sh 'translated_networks']);

%% test 10: load nn 5
cont5 = Load_nn([fp sh 'engine' sh 'nnmt' sh 'original_networks' sh 'ACASXU_run2a_5_8_batch_2000.nnet']);

%% test 11: NN TEST 6
cont6 = Load_nn([fp sh 'engine' sh 'nnmt' sh 'original_networks' sh 'CartPole_Controller.h5'],[fp sh 'engine' sh 'nnmt' sh 'translated_networks'],[fp sh 'engine' sh 'nnmt' sh 'original_networks' sh 'CartPole_Controller.json']);

%% test 12: load nn 7
cont7 = Load_nn([fp sh 'engine' sh 'nnmt' sh 'original_networks' sh 'example2' sh 'checkpoint'], 'models');
