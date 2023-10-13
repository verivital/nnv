%% Train all possible 2D classification models for medmnist data

medmnist_path = "data/mat_files/"; % path to data

datasets = dir(medmnist_path+"*.mat");

for i=1:length(datasets)
    if ~endsWith(datasets(i).name, "3d.mat") || ~contains(datasets(i).name, "chest")
        try
            train_medmnist2d(medmnist_path+datasets(i).name);
        catch ME
            warning("Failed!!")
            warning(ME.message);
            disp(medmnist_path+datasets(i).name)
        end
    end
end

% accuracy provided from the benchmarks provided (ResNets and more) can be
% found in https://medmnist.com/

% Best accuracies for each dataset
% --------------------------------------------------------
%    DATASET        BEST ACC        Small (ours)     Large (ours) 
% PathMNIST          0.911             0.8259
% ChestMNIST         0.948             0.9474
% DermaMNIST         0.768             0.7092
% OCTMNIST           0.776             0.6540
% PneumoniaMNIST     0.884             0.8189
% RetinaMNIST        0.531             0.4925
% BreastMNIST        0.863             0.8141
% BloodMNIST         0.966             0.8951
% TissueMNIST        0.681             0.5497
% OrganAMNIST        0.951             0.8613
% OrganCMNIST        0.920             0.8664
% OrganSMNIST        0.785             0.6952

