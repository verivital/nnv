%% Run examples to load a FNN from different formats using nnmt

sh = filesep;
% Get paths to the inputs for the first example
genp = split(string(pwd),sh);
rp = {'nnv','code','nnv','tests','io'};

if ~isequal(genp(end-4:end),string(rp)')
    error('Executing from %s. \nPlease, change your current folder path to ../nnv/code/nnv/io/tests and run script again',pwd);
end
% construct the oaths automatically
fp = genp(1:end-2);
fp = strjoin(fp,filesep);
fp = char(fp);
disp('Coverting neural network 1')
cont1 = Load_nn([fp sh 'engine' sh 'nnmt' sh 'original_networks' sh 'neural_network_information_13'],[fp sh 'engine' sh 'nnmt' sh 'translated_networks']);
disp('Coverting neural network 2')
cont2 = Load_nn([fp sh 'engine' sh 'nnmt' sh 'original_networks' sh 'neural_network_information_2'],[fp sh 'engine' sh 'nnmt' sh 'translated_networks']);
disp('Coverting neural network 3')
cont3 = Load_nn([fp sh 'engine' sh 'nnmt' sh 'original_networks' sh 'MC_16sigX16sigX1tanh.h5'],[fp sh 'engine' sh 'nnmt' sh 'translated_networks']);
disp('Coverting neural network 4')
cont4 = Load_nn([fp sh 'engine' sh 'nnmt' sh 'original_networks' sh 'ACASXU_run2a_4_3_batch_2000.nnet'],[fp sh 'engine' sh 'nnmt' sh 'translated_networks']);
disp('Coverting neural network 5')
cont5 = Load_nn([fp sh 'engine' sh 'nnmt' sh 'original_networks' sh 'ACASXU_run2a_5_8_batch_2000.nnet']);
disp('Coverting neural network 6')
cont6 = Load_nn([fp sh 'engine' sh 'nnmt' sh 'original_networks' sh 'CartPole_Controller.h5'],[fp sh 'engine' sh 'nnmt' sh 'translated_networks'],[fp sh 'engine' sh 'nnmt' sh 'original_networks' sh 'CartPole_Controller.json']);
disp('Coverting neural network 7')
cont7 = Load_nn([fp sh 'engine' sh 'nnmt' sh 'original_networks' sh 'example2' sh 'checkpoint'], 'models');

disp('We have succesfully converted all the neural networks!');
