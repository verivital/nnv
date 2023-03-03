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
acc_small = eval_small(Xall,Yall);
% Run medium network
acc_medium = eval_medium(Xall,Yall);
% Run larger network
acc_tiny = eval_tiny(Xall,Yall);
save('eval.mat','acc_small','acc_medium','acc_tiny');
