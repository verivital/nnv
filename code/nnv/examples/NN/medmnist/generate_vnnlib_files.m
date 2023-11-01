%% Let's create some examples for medmnist 2D and 3D

%% Dataset 1
% We start with bloodmnist
data1 = "data/mat_files/bloodmnist.mat";
load(data1);

test_images = permute(test_images, [2 3 4 1]);
test_labels = test_labels + 1;

outputSize = 8; % number of classes in dataset
N = 10; % number of vnnlib files to create
epsilon = [1,2,3]; % {epsilon} pixel color values for every channel

name = "vnnlib/bloodmnist_linf_";

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


%% Dataset 2
% Now we choose a grayscale dataset for 2D classification
% pneumoniamnist
data1 = "data/mat_files/pneumoniamnist.mat";
load(data1);

test_images = permute(test_images, [2 3 4 1]);
test_labels = test_labels + 1;

outputSize = 8; % number of classes in dataset
N = 10; % number of vnnlib files to create
epsilon = [1,2,3]; % {epsilon} pixel color values for every channel

name = "vnnlib/pneumoniamnist_linf_";

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
