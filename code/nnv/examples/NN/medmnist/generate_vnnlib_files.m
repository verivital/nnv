%% Let's create some examples for medmnist 2D and 3D

%% All 2D datasets

datapath = "data/mat_files/";
datafiles = ["bloodmnist"; "breastmnist.mat"; "dermamnist.mat"; "octmnist.mat"; "organamnist.mat"; ...
    "organcmnist.mat"; "organsmnist.mat"; "pathmnist.mat"; "pneumoniamnist"; "retinamnist.mat"; "tissuemnist.mat"];

N = 10; % number of vnnlib files to create
epsilon = [1,2,3]; % {epsilon} pixel color values for every channel

for k=1:length(datafiles)
    % load data
    load(datapath + datafiles(k));
    % preprocess dataa
    test_images = permute(test_images, [2 3 4 1]);
    test_labels = test_labels + 1;
    outputSize = length(unique(test_labels)); % number of classes in dataset
    % create file name
    dataname = split(datafiles(k), '.');
    name = "vnnlib/" + dataname{1} + "_linf_";
    % create vnnlib files
    for i=1:N
        img = test_images(:,:,:,i);
        outputSpec = create_output_spec(outputSize, test_labels(i));
        for j=1:length(epsilon)
            [lb,ub] = l_inf_attack(img, epsilon(j), 255, 0);
            vnnlibfile = name+string(epsilon(j))+"_"+string(i)+".vnnlib";
            export2vnnlib(lb, ub, outputSize, outputSpec, vnnlibfile);
            disp("Created property "+vnnlibfile);
        end
    end
end




%% Helper functions

% Return the bounds an linf attack
function [lb,ub] = l_inf_attack(img, epsilon, max_value, min_value)
    imgSize = size(img);
    disturbance = epsilon * ones(imgSize, "like", img); % disturbance value
    lb = max(img - disturbance, min_value);
    ub = min(img + disturbance, max_value);
    lb = single(lb);
    lb = reshape(lb, [], 1);
    ub = single(ub);
    ub = reshape(ub, [], 1);
end

% Define unsafe (not robust) property 
function Hs = create_output_spec(outSize, target)
    % @Hs: unsafe/not robust region defined as a HalfSpace
    %  - target: label idx of the given input set

    if target > outSize
        error("Target idx must be less than or equal to the output size of the NN.");
    end

    % Define HalfSpace Matrix and vector
    G = ones(outSize,1);
    G = diag(G);
    G(target, :) = [];
    G = -G;
    G(:, target) = 1;

    % Create HalfSapce to define robustness specification
    Hs = [];
    for i=1:height(G)
        Hs = [Hs; HalfSpace(G(i,:), 0)];
    end

end
