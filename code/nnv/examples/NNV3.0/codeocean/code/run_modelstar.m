function run_modelstar()
%% ModelStar Test Runner for CodeOcean
% Runs weight perturbation verification on MNIST MLP model
% Tests 3 layers (fc_6, fc_5, fc_4) to match original NNV 3.0 experiments
%
% Data required in /data/:
%   - MNIST/mnist_model_fc.mat
%
% Output in /results/ModelStar/:
%   - modelstar_results.csv
%   - modelstar_summary.txt

disp('Starting ModelStar (Weight Perturbation Verification)...');

%% Configuration for CodeOcean
config.modelDir = '/data/MNIST';
config.modelFile = 'mnist_model_fc.mat';
config.outputDir = '/results/ModelStar';
config.numImages = 100;         % Number of images to test (matches original)

% Layer configurations matching original MNIST_MLP.yaml
% Each layer has different perturbation fractions based on weight range
config.layers = {
    struct('name', 'fc_6', 'fracs', [0.005, 0.01, 0.015, 0.02]);
    struct('name', 'fc_5', 'fracs', [0.001, 0.002, 0.003, 0.004]);
    struct('name', 'fc_4', 'fracs', [0.001, 0.002, 0.003, 0.004]);
};

%% Validate paths
modelPath = fullfile(config.modelDir, config.modelFile);
if ~exist(modelPath, 'file')
    error('Model file not found: %s\nPlease upload mnist_model_fc.mat to /data/MNIST/', modelPath);
end

% Create output directory
if ~exist(config.outputDir, 'dir')
    mkdir(config.outputDir);
end

%% Load Network
disp('Loading MNIST MLP model...');
modelData = load(modelPath);
matlabnet = modelData.net;

% Create NNV model
net = matlab2nnv(matlabnet);
net.Name = 'MNIST_MLP';
net.OutputSize = 10;
netbs = getByteStreamFromArray(net);  % Byte stream for weight perturbation copies

disp('Model loaded successfully.');

%% Define Reachability Options
reachOptions = struct;
reachOptions.reachMethod = 'approx-star';
reachOptions.numCores = 1;
reachOptions.device = 'cpu';

%% Load MNIST Test Images
disp('Loading MNIST test images...');
[images, labels] = load_images_MNIST('database', 'mnist', 'n', config.numImages, 'matlabnet', matlabnet);
disp(['Loaded ' num2str(length(images)) ' images.']);

%% Initialize Results Storage
numLayers = length(config.layers);
allResults = cell(numLayers, 1);

%% Main Verification Loop - iterate over layers
disp(' ');
disp('========== Running ModelStar Verification ==========');
disp(['Testing ' num2str(numLayers) ' layers with ' num2str(config.numImages) ' images each']);

for layerNum = 1:numLayers
    layerConfig = config.layers{layerNum};
    layerName = layerConfig.name;
    pertFracs = layerConfig.fracs;

    disp(' ');
    disp(['====== Layer ' num2str(layerNum) '/' num2str(numLayers) ': ' layerName ' ======']);

    % Get layer index
    layerIdx = net.name2indx(layerName);
    disp(['Layer index: ' num2str(layerIdx) ', Perturbation fractions: ' mat2str(pertFracs)]);

    % Initialize results for this layer
    layerResults = struct();
    layerResults.name = layerName;
    layerResults.numImages = length(images);
    layerResults.fracs = pertFracs;
    layerResults.verified = zeros(1, length(pertFracs));
    layerResults.times = zeros(1, length(pertFracs));

    for fracIdx = 1:length(pertFracs)
        frac = pertFracs(fracIdx);
        disp(' ');
        disp(['--- Perturbation Fraction: ' num2str(frac) ' ---']);

        % Restore original network
        net = getArrayFromByteStream(netbs);

        % Get perturbation magnitude as fraction of weight range
        p = frac * WPutils.get_weights_range(net, layerIdx);

        % Apply perturbation to weights
        net.Layers{layerIdx}.perturb_whole_layer(-p, p);

        % Verify each image
        numVerified = 0;
        totalTime = 0;

        for imgIdx = 1:length(images)
            img = images{imgIdx};

            t = tic;
            try
                result = WPutils.verify_robustness_for_3dim_img(net, reachOptions, 'input', img);
                if result == 1
                    numVerified = numVerified + 1;
                    fprintf('+');
                else
                    fprintf('-');
                end
            catch ME
                fprintf('E');
                disp(['  Error on image ' num2str(imgIdx) ': ' ME.message]);
            end
            imgTime = toc(t);
            totalTime = totalTime + imgTime;
        end
        fprintf('\n');

        layerResults.verified(fracIdx) = numVerified;
        layerResults.times(fracIdx) = totalTime;

        percentVerified = 100 * numVerified / length(images);
        avgTime = totalTime / length(images);

        disp(['  Verified: ' num2str(numVerified) '/' num2str(length(images)) ' (' num2str(percentVerified, '%.1f') '%)']);
        disp(['  Total time: ' num2str(totalTime, '%.2f') 's, Avg per image: ' num2str(avgTime, '%.2f') 's']);
    end

    allResults{layerNum} = layerResults;
end

%% Save Results
disp(' ');
disp('========== Saving Results ==========');

% Save CSV with all layers
csvFile = fullfile(config.outputDir, 'modelstar_results.csv');
fid = fopen(csvFile, 'w');
fprintf(fid, 'Layer,PerturbationFrac,VerifiedCount,TotalImages,VerifiedPercent,TotalTime,AvgTimePerImage\n');

for layerNum = 1:numLayers
    layerResults = allResults{layerNum};
    for fracIdx = 1:length(layerResults.fracs)
        frac = layerResults.fracs(fracIdx);
        verified = layerResults.verified(fracIdx);
        totalTime = layerResults.times(fracIdx);
        percentVerified = 100 * verified / layerResults.numImages;
        avgTime = totalTime / layerResults.numImages;
        fprintf(fid, '%s,%.6f,%d,%d,%.2f,%.4f,%.4f\n', ...
            layerResults.name, frac, verified, layerResults.numImages, percentVerified, totalTime, avgTime);
    end
end
fclose(fid);
disp(['Saved: ' csvFile]);

% Save Summary
summaryFile = fullfile(config.outputDir, 'modelstar_summary.txt');
fid = fopen(summaryFile, 'w');
fprintf(fid, 'ModelStar Weight Perturbation Verification Results\n');
fprintf(fid, '==================================================\n\n');
fprintf(fid, 'Model: MNIST_MLP\n');
fprintf(fid, 'Number of layers tested: %d\n', numLayers);
fprintf(fid, 'Number of images per test: %d\n\n', config.numImages);

for layerNum = 1:numLayers
    layerResults = allResults{layerNum};
    fprintf(fid, 'Layer: %s\n', layerResults.name);
    fprintf(fid, '--------\n');
    for fracIdx = 1:length(layerResults.fracs)
        frac = layerResults.fracs(fracIdx);
        verified = layerResults.verified(fracIdx);
        percentVerified = 100 * verified / layerResults.numImages;
        fprintf(fid, '  Perturbation %.4f: %d/%d (%.1f%%) verified\n', ...
            frac, verified, layerResults.numImages, percentVerified);
    end
    fprintf(fid, '\n');
end
fprintf(fid, 'Date: %s\n', datestr(now));
fclose(fid);
disp(['Saved: ' summaryFile]);

%% Display Summary
disp(' ');
disp('========== ModelStar Summary ==========');
disp(['Model: MNIST_MLP']);
disp(['Layers tested: ' num2str(numLayers)]);
disp(['Images per test: ' num2str(config.numImages)]);
disp(' ');

for layerNum = 1:numLayers
    layerResults = allResults{layerNum};
    disp(['Layer: ' layerResults.name]);
    disp('Perturbation  | Verified  | Percentage');
    disp('--------------|-----------|------------');
    for fracIdx = 1:length(layerResults.fracs)
        frac = layerResults.fracs(fracIdx);
        verified = layerResults.verified(fracIdx);
        percentVerified = 100 * verified / layerResults.numImages;
        disp(sprintf('   %.4f     |    %2d     |   %.1f%%', frac, verified, percentVerified));
    end
    disp(' ');
end

disp('ModelStar verification complete.');

end
