% This script is used for generating .vnnlib files for the UCF11 verification experiments described in the SoSym journal extension
% of "Robustness Verification for Video Classification Neural Networks". 
%
% This script is to be used with defined epsilons, number of frames, and the samples to generate the vnnlib files for.
% i.e., if you want to generate vnnlib files for both the 8-frame and 16-frame versions of the UCF11 verification dataset,
% then you will need to run this script twice.

% set hyperparameters
epsilons = [1/255, 2/255, 3/255];
numFrames = 32;

% load the data
data = readNPY(sprintf("data/UCF11/ucf11_grayscale_%df_verification_data.npy", numFrames)); % 100 numFrames 112 112 1
labels = readNPY(sprintf("data/UCF11/ucf11_grayscale_%df_verification_labels.npy", numFrames));

numSamples = 9; % number of samples whose vnnlib files we need to generate
% in this case, we do 1-9 because this is all we have verified for UCF11 so far...
% assuming we did more, we would want to go from 10 -> N where N is the next total we choose.

% path to vnnlib directory
vnnlibpath = "/tmp/vnnlib/";

% loop
for sampleIndex=1:numSamples
    fprintf("Iteration: %d\n", sampleIndex)
    for epsIndex=1:length(epsilons)
        s = data(sampleIndex,:,:,:);
        s = squeeze(s); % need to squeeze to match dimensions
        l = labels(sampleIndex)+1;
        e = epsilons(epsIndex);
        generate_vnnlib(s, l, sampleIndex, e, epsIndex, numFrames, vnnlibpath);
    end

end

fprintf("Done generating .vnnlib files for samples from UCF11 with %d frames!", numFrames);

%% Helper Functions
function generate_vnnlib(sample, label, index, epsilon, epsIndex, numFrames, vnnlibpath)
    [VS, lb, ub] = L_inf_attack(sample, epsilon, numFrames);
    outputSize = 11; % there are 11 classes
    outputSpec = create_output_spec(outputSize, label);
    vnnlibfile = vnnlibpath + sprintf("ucf11/%d/%d_255/sample_%d.vnnlib", numFrames, epsIndex, index);
    export2vnnlib(lb, ub, outputSize, outputSpec, vnnlibfile);
end

function [VS, lb_clip, ub_clip] = L_inf_attack(x, epsilon, numFrames)
    lb = squeeze(x);
    ub = squeeze(x);

    % Perturb the frames
    for fn=1:numFrames
        lb(fn, :, :, :) = x(fn, :, :, :) - epsilon;
        ub(fn, :, :, :) = x(fn, :, :, :) + epsilon;
    end

    % Clip the perturbed values to be between 0-1
    lb_min = zeros(numFrames, 112, 112, 1);
    ub_max = ones(numFrames, 112, 112, 1);
    lb_clip = max(lb, lb_min);
    ub_clip = min(ub, ub_max);

    % Create the volume star
    VS = VolumeStar(lb_clip, ub_clip);
end

function Hs = create_output_spec(outSize, target)
    if target > outSize
        error("Target idx must be less than or equal to the output size of the NN")
    end

    % Define HalfSpace matrix and vector
    G = ones(outSize, 1);
    G = diag(G);
    G(target, :) = [];
    G = -G;
    G(:, target) = 1;

    % Create HalfSpace to define robustness specification
    Hs = [];
    for i=1:height(G)
        Hs = [Hs; HalfSpace(G(i,:), 0)];
    end
end
