%% Verify all possible 2D classification models for medmnist data

medmnist_path = "data/mat_files/"; % path to data

datasets = dir(medmnist_path+"*.mat");

for i=1:length(datasets)
    if ~endsWith(datasets(i).name, "3d.mat")
        try
            verify_medmnist(medmnist_path+datasets(i).name);
        catch ME
            warning("Failed!!")
            warning(ME.message);
            disp(medmnist_path+datasets(i).name)
        end
    end
end