%% run_probver.m - YOLO Benchmark Verification Script
% Standalone script for running probabilistic verification on the yolo_2023 benchmark
% using the cp-star reachability method.
%
% PREREQUISITES: Run startup_nnv.m from the NNV root directory first, or
% this script will attempt to do it automatically.

%% ========== CONFIGURATION ==========
numSamples = 3;   % Number of instances to verify
randomSeed = 42;  % Seed for reproducibility
nRand = 100;      % Number of random samples for falsification

%% ========== SETUP PATHS ==========
scriptDir = fileparts(mfilename('fullpath'));
benchmarkDir = fullfile(scriptDir, 'yolo_2023');
onnxDir = fullfile(benchmarkDir, 'onnx');
vnnlibDir = fullfile(benchmarkDir, 'vnnlib');
instancesFile = fullfile(benchmarkDir, 'instances.csv');
resultsFile = fullfile(scriptDir, 'results_summary.csv');

%% ========== DECOMPRESS FILES IF NEEDED ==========
disp('Checking for compressed files...');

% Decompress ONNX model if needed
onnxFile = fullfile(onnxDir, 'TinyYOLO.onnx');
onnxGzFile = fullfile(onnxDir, 'TinyYOLO.onnx.gz');
if ~isfile(onnxFile) && isfile(onnxGzFile)
    disp('Decompressing TinyYOLO.onnx.gz...');
    gunzip(onnxGzFile, onnxDir);
end

% Decompress all vnnlib files if needed
vnnlibGzFiles = dir(fullfile(vnnlibDir, '*.vnnlib.gz'));
for i = 1:length(vnnlibGzFiles)
    gzFile = fullfile(vnnlibDir, vnnlibGzFiles(i).name);
    uncompressedFile = fullfile(vnnlibDir, erase(vnnlibGzFiles(i).name, '.gz'));
    if ~isfile(uncompressedFile)
        disp(['Decompressing ' vnnlibGzFiles(i).name '...']);
        gunzip(gzFile, vnnlibDir);
    end
end

%% ========== PARSE INSTANCES.CSV ==========
disp('Reading instances.csv...');
fid = fopen(instancesFile, 'r');
instances = {};
while ~feof(fid)
    line = fgetl(fid);
    if ischar(line) && ~isempty(strtrim(line))
        parts = strsplit(line, ',');
        if length(parts) >= 3
            instances{end+1} = struct(...
                'onnx', strtrim(parts{1}), ...
                'vnnlib', strtrim(parts{2}), ...
                'timeout', str2double(parts{3}));
        end
    end
end
fclose(fid);

totalInstances = length(instances);
disp(['Found ' num2str(totalInstances) ' instances']);

% Randomly select instances with seeded RNG for reproducibility
rng(randomSeed);
numSamples = min(numSamples, totalInstances);
selectedIndices = randperm(totalInstances, numSamples);
selectedIndices = sort(selectedIndices);  % Sort for cleaner output
disp(['Randomly selected instances (seed=' num2str(randomSeed) '): ' mat2str(selectedIndices)]);

%% ========== LOAD NETWORK (once for all instances) ==========
disp('Loading TinyYOLO network...');
onnxPath = fullfile(benchmarkDir, 'onnx', 'TinyYOLO.onnx');
net = importNetworkFromONNX(onnxPath, "InputDataFormats", "BCSS");

try
    nnvnet = matlab2nnv(net);
catch
    nnvnet = "";
end

inputSize = net.Layers(1, 1).InputSize;
inputFormat = "default";
needReshape = 2;

% YOLO-specific reachability options (cp-star)
reachOptions = struct;
reachOptions.train_epochs = 200;
reachOptions.train_lr = 0.0001;
reachOptions.coverage = 0.999;
reachOptions.confidence = 0.999;
reachOptions.train_mode = 'Linear';
reachOptions.surrogate_dim = [-1, -1];
reachOptions.threshold_normal = 1e-5;
reachOptions.reachMethod = 'cp-star';
reachOptionsList{1} = reachOptions;

disp('Network loaded successfully.');
disp(['Input size: ' mat2str(inputSize)]);

%% ========== MAIN VERIFICATION LOOP ==========
results = cell(numSamples, 1);
totalStartTime = tic;

for i = 1:numSamples
    idx = selectedIndices(i);
    instance = instances{idx};

    disp(' ');
    disp(['========== Instance ' num2str(idx) '/' num2str(totalInstances) ' (' num2str(i) '/' num2str(numSamples) ') ==========']);
    disp(['VNNLIB: ' instance.vnnlib]);

    result = struct();
    result.index = idx;
    result.onnx = instance.onnx;
    result.vnnlib = instance.vnnlib;
    result.status = 'unknown';
    result.statusCode = 2;  % 0=sat, 1=unsat, 2=unknown, -1=error
    result.time = 0;
    result.error = '';

    try
        t = tic;
        status = 2; % unknown

        % Load property to verify
        vnnlibPath = fullfile(benchmarkDir, instance.vnnlib);
        property = load_vnnlib(vnnlibPath);
        lb = property.lb;
        ub = property.ub;
        prop = property.prop;

        %% Falsification (SAT check)
        counterEx = nan;
        if ~isa(lb, "cell") && length(prop) == 1
            counterEx = falsify_single(net, lb, ub, inputSize, nRand, prop{1}.Hg, needReshape, inputFormat);
        elseif isa(lb, "cell") && length(lb) == length(prop)
            for spc = 1:length(lb)
                counterEx = falsify_single(net, lb{spc}, ub{spc}, inputSize, nRand, prop{spc}.Hg, needReshape, inputFormat);
                if iscell(counterEx)
                    break
                end
            end
        elseif isa(lb, "cell") && length(prop) == 1
            for arr = 1:length(lb)
                counterEx = falsify_single(net, lb{arr}, ub{arr}, inputSize, nRand, prop{1}.Hg, needReshape, inputFormat);
                if iscell(counterEx)
                    break
                end
            end
        end

        cEX_time = toc(t);

        if iscell(counterEx)
            status = 0; % SAT (property violated)
            disp(['Counterexample found in ' num2str(cEX_time) 's']);
        else
            %% Reachability (UNSAT check)
            vT = tic;

            if ~isa(lb, "cell") && length(prop) == 1
                if ~nnz(lb-ub)
                    status = 1; % verified (point input)
                else
                    IS = create_input_set(lb, ub, inputSize, needReshape);
                    ySet = Prob_reach(net, IS, reachOptions);
                    status = verify_specification(ySet, prop);
                end
            elseif isa(lb, "cell") && length(lb) == length(prop)
                local_status = 2*ones(length(lb), 1);
                for spc = 1:length(lb)
                    lb_spc = lb{spc};
                    ub_spc = ub{spc};
                    if ~nnz(lb_spc - ub_spc)
                        local_status(spc) = 1;
                    else
                        IS = create_input_set(lb_spc, ub_spc, inputSize, needReshape);
                        ySet = Prob_reach(net, IS, reachOptions);
                        if isempty(ySet.C)
                            dd = ySet.V; DD = ySet.V;
                            ySet = Star(dd, DD);
                        end
                        local_status(spc) = verify_specification(ySet, prop(spc));
                    end
                end
                if all(local_status == 1)
                    status = 1;
                else
                    status = 2;
                end
            elseif isa(lb, "cell") && length(prop) == 1
                local_status = 2*ones(length(lb), 1);
                for spc = 1:length(lb)
                    lb_spc = lb{spc};
                    ub_spc = ub{spc};
                    if ~nnz(lb_spc - ub_spc)
                        local_status(spc) = 1;
                    else
                        IS = create_input_set(lb_spc, ub_spc, inputSize, needReshape);
                        ySet = Prob_reach(net, IS, reachOptions);
                        tempStatus = verify_specification(ySet, prop(spc));
                        local_status(spc) = tempStatus;
                    end
                end
                if all(local_status == 1)
                    status = 1;
                else
                    status = 2;
                end
            end

            vT = toc(vT);
            disp(['Reachability time: ' num2str(vT) 's']);
        end

        tTime = toc(t);

        % Update result
        result.statusCode = status;
        if status == 0
            result.status = 'sat';
        elseif status == 1
            result.status = 'unsat';
        else
            result.status = 'unknown';
        end
        result.time = tTime;

        disp(['Result: ' result.status ' (' num2str(tTime) 's)']);

    catch ME
        result.statusCode = -1;
        result.status = 'error';
        result.error = ME.message;
        result.time = toc(t);
        disp(['ERROR: ' ME.message]);
    end

    results{i} = result;
end

totalTime = toc(totalStartTime);

%% ========== SAVE AND DISPLAY RESULTS ==========
disp(' ');
disp('========== RESULTS SUMMARY ==========');

% Count results
satCount = 0;
unsatCount = 0;
unknownCount = 0;
errorCount = 0;

for i = 1:length(results)
    r = results{i};
    switch r.status
        case 'sat'
            satCount = satCount + 1;
        case 'unsat'
            unsatCount = unsatCount + 1;
        case 'unknown'
            unknownCount = unknownCount + 1;
        case 'error'
            errorCount = errorCount + 1;
    end
end

disp(['Total instances: ' num2str(length(results))]);
disp(['  SAT (violated):   ' num2str(satCount)]);
disp(['  UNSAT (verified): ' num2str(unsatCount)]);
disp(['  Unknown:          ' num2str(unknownCount)]);
disp(['  Errors:           ' num2str(errorCount)]);
disp(['Total time: ' num2str(totalTime) 's']);

% Save to CSV
disp(' ');
disp(['Saving results to ' resultsFile '...']);
fid = fopen(resultsFile, 'w');
fprintf(fid, 'index,onnx,vnnlib,status,time,error\n');
for i = 1:length(results)
    r = results{i};
    % Escape any commas in error message
    errorMsg = strrep(r.error, ',', ';');
    errorMsg = strrep(errorMsg, newline, ' ');
    fprintf(fid, '%d,%s,%s,%s,%.4f,%s\n', r.index, r.onnx, r.vnnlib, r.status, r.time, errorMsg);
end
fclose(fid);
disp('Results saved.');

% Display errors if any
if errorCount > 0
    disp(' ');
    disp('========== ERRORS ==========');
    for i = 1:length(results)
        r = results{i};
        if strcmp(r.status, 'error')
            disp(['Instance ' num2str(r.index) ': ' r.error]);
        end
    end
end

%% ========== HELPER FUNCTIONS ==========

function IS = create_input_set(lb, ub, inputSize, needReshape)
    % Get input bounds
    if ~isscalar(inputSize)
        lb = reshape(lb, inputSize);
        ub = reshape(ub, inputSize);
    end

    % Format bounds into correct dimensions
    if needReshape == 1
        lb = permute(lb, [2 1 3]);
        ub = permute(ub, [2 1 3]);
    elseif needReshape == 2
        newSize = [inputSize(2) inputSize(1) inputSize(3:end)];
        lb = reshape(lb, newSize);
        lb = permute(lb, [2 1 3 4]);
        ub = reshape(ub, newSize);
        ub = permute(ub, [2 1 3 4]);
    end

    % Create input set
    IS = ImageStar(lb, ub);

    % Delete constraints for non-interval dimensions
    try
        xxx = find((lb-ub));
        xxx = reshape(xxx, [], 1);
        if numel(lb) ~= length(xxx)
            IS.C = IS.C(:,xxx);
            IS.pred_lb = IS.pred_lb(xxx);
            IS.pred_ub = IS.pred_ub(xxx);
            xxx = xxx + 1;
            IS.V = IS.V(:,:,:,[1;xxx]);
            IS.numPred = length(xxx);
        end
    end
end

function xRand = create_random_examples(net, lb, ub, nR, inputSize, needReshape, inputFormat)
    xB = Box(lb, ub);
    xRand = xB.sample(nR-2);
    xRand = [lb, ub, xRand];
    if needReshape
        if needReshape == 2
            newSize = [inputSize(2) inputSize(1) inputSize(3:end)];
            xRand = reshape(xRand, [newSize nR]);
            xRand = permute(xRand, [2 1 3 4]);
        else
            xRand = reshape(xRand, [inputSize nR]);
            xRand = permute(xRand, [2 1 3 4]);
        end
    else
        xRand = reshape(xRand, [inputSize nR]);
    end
    if isa(net, 'dlnetwork')
        if strcmp(inputFormat, "default")
            if isa(net.Layers(1, 1), 'nnet.cnn.layer.ImageInputLayer')
                xRand = dlarray(xRand, "SSCB");
            elseif isa(net.Layers(1, 1), 'nnet.cnn.layer.FeatureInputLayer') || isa(net.Layers(1, 1), 'nnet.onnx.layer.FeatureInputLayer')
                xRand = dlarray(xRand, "CB");
            else
                disp(net.Layers(1,1));
                error("Unknown input format");
            end
        else
            if contains(inputFormat, "U")
                xRand = dlarray(xRand, inputFormat+"U");
            else
                xRand = dlarray(xRand, inputFormat);
            end
        end
    end
end

function counterEx = falsify_single(net, lb, ub, inputSize, nRand, Hs, needReshape, inputFormat)
    counterEx = nan;
    xRand = create_random_examples(net, lb, ub, nRand, inputSize, needReshape, inputFormat);
    s = size(xRand);
    n = length(s);
    for i = 1:s(n)
        x = get_example(xRand, i);
        try
            yPred = predict(net, x);
            if isa(x, 'dlarray')
                x = extractdata(x);
                yPred = extractdata(yPred);
            end
            yPred = reshape(yPred, [], 1);
            for h = 1:length(Hs)
                if Hs(h).contains(double(yPred))
                    n = numel(x);
                    if needReshape == 2
                        x = permute(x, [2 1 3]);
                    elseif needReshape == 1
                        if ndims(x) == 3
                            x = permute(x, [2 1 3]);
                        end
                    end
                    counterEx = {x; yPred};
                    return;
                end
            end
        end
    end
end

function x = get_example(xRand, i)
    s = size(xRand);
    n = length(s);
    if n == 4
        x = xRand(:,:,:,i);
    elseif n == 3
        x = xRand(:,:,i);
    elseif n == 2
        x = xRand(:,i);
        xsize = size(x);
        if xsize(1) ~= 1 && ~isa(x, "dlarray")
            x = x';
        end
    else
        error("InputSize = " + string(s));
    end
end
