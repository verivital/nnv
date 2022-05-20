% Test all CNNs
% Load all test images
Xall = processMNISTimages('t10k-images.idx3-ubyte');
Yall = processMNISTlabels('t10k-labels.idx1-ubyte');
Xall = extractdata(Xall);
Yall = double(Yall);
% Use full test set
numT = length(Yall);
% Use 100 to test the test functions
% Xall = Xall(:,:,:,1:100);
% Yall = Yall(1:100,:);
% Run smaller network
acc_small = eval_ffnn_small(Xall,Yall);
% Run medium network
acc_medium = eval_ffnn(Xall,Yall);
% Run larger network
acc_large = eval_ffnn_large(Xall,Yall);
save('eval.mat','acc_small','acc_medium','acc_large');
