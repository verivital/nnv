%% Train all possible 3D classification models for medmnist data

medmnist_path = "data/mat_files/"; % path to data

datasets = dir(medmnist_path+"*.mat");

for i=1:length(datasets)
    if endsWith(datasets(i).name, "3d.mat")
        try
            train_medmnist3d(medmnist_path+datasets(i).name);
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
% ------------------------------------------------
%    DATASET           BEST ACC        Small (ours)   
% AdrenalMNIST3d        0.802             0.788
% FractureMNIST3d       0.517             0.450
% NoduleMNIST3d         0.848             0.848
% OrganMNIST3d          0.907             0.874
% SynapseMNIST3d        0.745             0.730
% VesselMNIST3d         0.928             0.903


