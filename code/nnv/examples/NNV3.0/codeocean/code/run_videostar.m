function run_videostar()
%% VideoStar Test Runner for CodeOcean
% Runs video classification verification using VideoStar
%
% Data required in /data/:
%   - FORMALISE2025/models_mat/zoomin_4f.mat (pre-converted network)
%   - FORMALISE2025/data/ZoomIn/*.npy
%
% Output in /results/VideoStar/:
%   - eps=1_255.csv, eps=2_255.csv, eps=3_255.csv

disp('Starting VideoStar (Video Classification Verification)...');

%% Suppress warnings to clean up output
warning('off', 'all');

%% Configuration for CodeOcean
config.matDir = '/data/FORMALISE2025/models_mat';
config.dataDir = '/data/FORMALISE2025/data/ZoomIn';
config.outputDir = '/results/VideoStar';
config.dsType = 'zoom_in';          % Dataset type
config.sampleLen = 4;               % Number of frames
config.verAlgorithm = 'relax';      % Verification algorithm
config.numClasses = 10;
config.epsilon = [1/255; 2/255; 3/255];
config.timeout = 1800;              % 30 minutes per sample
config.sampleIndices = 1:10;        % Subset for quick testing

%% Validate paths
if ~exist(config.matDir, 'dir')
    error('Models directory not found: %s\nPlease upload FORMALISE2025/models_mat to /data/', config.matDir);
end

if ~exist(config.dataDir, 'dir')
    error('Data directory not found: %s\nPlease upload FORMALISE2025/data/ZoomIn to /data/', config.dataDir);
end

% Check if model exists
modelName = sprintf('zoomin_%df.mat', config.sampleLen);
modelPath = fullfile(config.matDir, modelName);
if ~exist(modelPath, 'file')
    error('Model not found: %s\nPlease create using convert_onnx_to_mat.m', modelPath);
end
disp(['Found model: ' modelName]);

% Create output directory
if ~exist(config.outputDir, 'dir')
    mkdir(config.outputDir);
end

%% Load Data
disp('Loading video data...');
dataFile = fullfile(config.dataDir, sprintf('mnistvideo_%s_%df_test_data_seq.npy', config.dsType, config.sampleLen));
labelsFile = fullfile(config.dataDir, sprintf('mnistvideo_%s_test_labels_seq.npy', config.dsType));

if ~exist(dataFile, 'file')
    error('Data file not found: %s', dataFile);
end
if ~exist(labelsFile, 'file')
    error('Labels file not found: %s', labelsFile);
end

data = readNPY(dataFile);
labels = readNPY(labelsFile);
disp(['Loaded data: ' num2str(size(data, 1)) ' samples']);

% Preprocessing
reshaped_data = permute(data, [1, 3, 2, 4, 5]); % to match BCSSS
data_squeezed = squeeze(reshaped_data);

%% Load Model from pre-converted .mat file
disp('Loading model...');
loaded = load(modelPath);
net = loaded.net;  % NNV network
net.OutputSize = config.numClasses;
disp(['Loaded pre-converted network: ' modelName]);

%% Verification Settings
reachOptions = struct;
if strcmp(config.verAlgorithm, 'relax')
    reachOptions.reachMethod = 'relax-star-area';
    reachOptions.relaxFactor = 0.5;
elseif strcmp(config.verAlgorithm, 'approx')
    reachOptions.reachMethod = 'approx-star';
end

%% Initialize Results Storage
numSamples = length(config.sampleIndices);
numEpsilons = length(config.epsilon);
results = cell(numSamples, numEpsilons);
times = cell(numSamples, numEpsilons);
methods = cell(numSamples, numEpsilons);

%% Run Verification
disp(' ');
disp('========== Running VideoStar Verification ==========');
disp(['Dataset: ' config.dsType]);
disp(['Sample Length: ' num2str(config.sampleLen) ' frames']);
disp(['Algorithm: ' config.verAlgorithm]);
disp(['Samples: ' num2str(numSamples)]);
disp('Epsilon values: 1/255, 2/255, 3/255');
disp(' ');

for sampleIdx = 1:numSamples
    sampleNum = config.sampleIndices(sampleIdx);
    disp(['Processing sample ' num2str(sampleNum) ' (' num2str(sampleIdx) '/' num2str(numSamples) ')...']);

    for epsIdx = 1:numEpsilons
        eps = config.epsilon(epsIdx);
        disp(['  Epsilon = ' num2str(epsIdx) '/255...']);

        try
            % Get the sample
            sample = squeeze(data_squeezed(sampleNum,:,:,:));
            label = labels(sampleNum) + 1;  % MATLAB is 1-indexed

            % Perform L_inf attack
            VS = L_inf_attack(sample, eps, config.sampleLen);

            t = tic;
            met = config.verAlgorithm;

            % Run verification
            res = net.verify_robustness(VS, reachOptions, label);

            results{sampleIdx, epsIdx} = res;
            times{sampleIdx, epsIdx} = toc(t);
            methods{sampleIdx, epsIdx} = met;

            disp(['    Result: ' num2str(res) ', Time: ' num2str(times{sampleIdx, epsIdx}, '%.2f') 's']);

        catch ME
            disp(['    Error: ' ME.message]);
            results{sampleIdx, epsIdx} = -1;
            times{sampleIdx, epsIdx} = -1;
            methods{sampleIdx, epsIdx} = ME.message;
        end
    end
end

%% Save Results
disp(' ');
disp('========== Saving Results ==========');

% Save results for each epsilon value
for epsIdx = 1:numEpsilons
    outputFile = fullfile(config.outputDir, sprintf('eps=%d_255.csv', epsIdx));

    fid = fopen(outputFile, 'w');
    fprintf(fid, 'SampleNumber,Result,Time,Method\n');

    for sampleIdx = 1:numSamples
        sampleNum = config.sampleIndices(sampleIdx);
        res = results{sampleIdx, epsIdx};
        t = times{sampleIdx, epsIdx};
        met = methods{sampleIdx, epsIdx};

        if isnumeric(t)
            fprintf(fid, '%d,%d,%.6f,%s\n', sampleNum, res, t, met);
        else
            fprintf(fid, '%d,%d,%s,%s\n', sampleNum, res, t, met);
        end
    end

    fclose(fid);
    disp(['Saved: ' outputFile]);
end

%% Display Summary
disp(' ');
disp('========== VideoStar Summary ==========');

for epsIdx = 1:numEpsilons
    verified = 0;
    unknown = 0;
    violated = 0;
    errors = 0;
    totalTime = 0;
    validTimes = 0;

    for sampleIdx = 1:numSamples
        res = results{sampleIdx, epsIdx};
        t = times{sampleIdx, epsIdx};

        if res == 1
            verified = verified + 1;
        elseif res == 0
            violated = violated + 1;
        elseif res == 2
            unknown = unknown + 1;
        else
            errors = errors + 1;
        end

        if isnumeric(t) && t > 0
            totalTime = totalTime + t;
            validTimes = validTimes + 1;
        end
    end

    avgTime = totalTime / max(validTimes, 1);

    disp(['Epsilon = ' num2str(epsIdx) '/255:']);
    disp(['  Verified (robust): ' num2str(verified)]);
    disp(['  Violated: ' num2str(violated)]);
    disp(['  Unknown: ' num2str(unknown)]);
    disp(['  Errors: ' num2str(errors)]);
    disp(['  Average time: ' num2str(avgTime, '%.2f') 's']);
    disp(' ');
end

disp('VideoStar verification complete.');
disp(['Results saved to: ' config.outputDir]);

end

%% Helper Function: L_inf attack
function VS = L_inf_attack(x, epsilon, numFrames)
    lb = squeeze(x);
    ub = squeeze(x);

    % Perturb the frames
    for fn=1:numFrames
        lb(fn, :, :) = x(fn, :, :) - epsilon;
        ub(fn, :, :) = x(fn, :, :) + epsilon;
    end

    % Clip the perturbed values to be between 0-1
    lb_min = zeros(numFrames, 32, 32);
    ub_max = ones(numFrames, 32, 32);
    lb_clip = max(lb, lb_min);
    ub_clip = min(ub, ub_max);

    % Create the volume star
    VS = VolumeStar(lb_clip, ub_clip);
end
