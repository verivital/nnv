function run_probver()
%% ProbVer Test Runner for CodeOcean
% Runs probabilistic verification on YOLO benchmark using CP-Star method
% Based on original: /code/nnv/examples/NNV3.0/ProbVer/run_probver.m
%
% Data required in /data/:
%   - ProbVer/yolo_2023/onnx/TinyYOLO.onnx (ONNX model)
%   - ProbVer/yolo_2023/vnnlib/*.vnnlib.gz
%   - ProbVer/yolo_2023/instances.csv
%
% Output in /results/ProbVer/:
%   - results_summary.csv

disp('Starting ProbVer (Probabilistic Verification)...');

%% Configuration
config.benchmarkDir = '/data/ProbVer/yolo_2023';
config.outputDir = '/results/ProbVer';
config.numSamples = 3;      % Number of instances to verify
config.randomSeed = 42;     % Seed for reproducibility
config.nRand = 100;         % Number of random samples for falsification

%% Validate paths
if ~exist(config.benchmarkDir, 'dir')
    error('Benchmark directory not found: %s', config.benchmarkDir);
end

onnxDir = fullfile(config.benchmarkDir, 'onnx');
vnnlibDir = fullfile(config.benchmarkDir, 'vnnlib');
instancesFile = fullfile(config.benchmarkDir, 'instances.csv');

if ~exist(onnxDir, 'dir')
    error('ONNX directory not found: %s', onnxDir);
end
if ~exist(vnnlibDir, 'dir')
    error('VNNLIB directory not found: %s', vnnlibDir);
end
if ~exist(instancesFile, 'file')
    error('instances.csv not found: %s', instancesFile);
end

if ~exist(config.outputDir, 'dir')
    mkdir(config.outputDir);
end

resultsFile = fullfile(config.outputDir, 'results_summary.csv');

%% Decompress VNNLIB files if needed
vnnlibGzFiles = dir(fullfile(vnnlibDir, '*.vnnlib.gz'));
for i = 1:length(vnnlibGzFiles)
    gzFile = fullfile(vnnlibDir, vnnlibGzFiles(i).name);
    uncompressedFile = fullfile(vnnlibDir, erase(vnnlibGzFiles(i).name, '.gz'));
    if ~isfile(uncompressedFile)
        gunzip(gzFile, vnnlibDir);
    end
end

%% Parse instances.csv
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
rng(config.randomSeed);
numSamples = min(config.numSamples, totalInstances);
selectedIndices = sort(randperm(totalInstances, numSamples));
disp(['Testing ' num2str(numSamples) ' of ' num2str(totalInstances) ' instances']);

%% Load network
disp('Loading TinyYOLO network...');
onnxPath = fullfile(onnxDir, 'TinyYOLO.onnx');

if ~exist(onnxPath, 'file')
    onnxGzPath = [onnxPath '.gz'];
    if exist(onnxGzPath, 'file')
        gunzip(onnxGzPath, onnxDir);
    else
        error('ONNX file not found: %s', onnxPath);
    end
end

try
    net = importNetworkFromONNX(onnxPath, "InputDataFormats", "BCSS");
    inputSize = net.Layers(1).InputSize;

    % Remove ONNX-specific layers that cause parallel pool serialization issues
    try
        lg = layerGraph(net);
        layersToRemove = {};
        for i = 1:numel(net.Layers)
            layerClass = class(net.Layers(i));
            if contains(layerClass, 'nnet.onnx.layer')
                layersToRemove{end+1} = net.Layers(i).Name;
            end
        end

        for i = 1:length(layersToRemove)
            layerName = layersToRemove{i};
            conns = lg.Connections;
            srcLayers = conns.Source(strcmp(conns.Destination, layerName));
            dstLayers = conns.Destination(strcmp(conns.Source, layerName));
            lg = removeLayers(lg, layerName);
            if ~isempty(srcLayers) && ~isempty(dstLayers)
                for j = 1:length(dstLayers)
                    try
                        lg = connectLayers(lg, srcLayers{1}, dstLayers{j});
                    catch
                    end
                end
            end
        end
        net = dlnetwork(lg);
        inputSize = net.Layers(1).InputSize;
    catch
        % Proceed with original network
    end

catch ME
    % Try fallback
    try
        net = importONNXNetwork(onnxPath);
        inputSize = net.Layers(1).InputSize;
    catch ME2
        error('Could not load ONNX network: %s', ME.message);
    end
end

disp(['Network loaded: ' num2str(numel(net.Layers)) ' layers, input size ' mat2str(inputSize)]);

inputFormat = "default";
needReshape = 2;

%% Reachability options (CP-Star)
reachOptions = struct;
reachOptions.train_epochs = 200;
reachOptions.train_lr = 0.0001;
reachOptions.coverage = 0.999;
reachOptions.confidence = 0.999;
reachOptions.train_mode = 'Linear';
reachOptions.surrogate_dim = [-1, -1];
reachOptions.threshold_normal = 1e-5;
reachOptions.reachMethod = 'cp-star';
reachOptions.train_device = 'cpu';
reachOptions.numCores = 1;  % Disable parallel to avoid ONNX layer serialization issues

%% Disable parallel pool (ONNX networks cannot be serialized to workers)
try
    poolobj = gcp('nocreate');
    if ~isempty(poolobj)
        delete(poolobj);
    end
    ps = parallel.Settings;
    ps.Pool.AutoCreate = false;
catch
end

%% Main verification loop
disp(' ');
disp('Running verification...');

results = cell(numSamples, 1);
totalStartTime = tic;

for i = 1:numSamples
    idx = selectedIndices(i);
    instance = instances{idx};

    result = struct();
    result.index = idx;
    result.onnx = instance.onnx;
    result.vnnlib = instance.vnnlib;
    result.status = 'unknown';
    result.statusCode = 2;
    result.time = 0;
    result.error = '';

    try
        t = tic;
        status = 2;

        % Load property
        vnnlibPath = fullfile(config.benchmarkDir, instance.vnnlib);
        if ~exist(vnnlibPath, 'file')
            error('VNNLIB file not found: %s', vnnlibPath);
        end

        property = load_vnnlib(vnnlibPath);
        lb = property.lb;
        ub = property.ub;
        prop = property.prop;

        % Falsification
        counterEx = nan;
        if ~isa(lb, "cell") && length(prop) == 1
            counterEx = falsify_single(net, lb, ub, inputSize, config.nRand, prop{1}.Hg, needReshape, inputFormat);
        elseif isa(lb, "cell") && length(lb) == length(prop)
            for spc = 1:length(lb)
                counterEx = falsify_single(net, lb{spc}, ub{spc}, inputSize, config.nRand, prop{spc}.Hg, needReshape, inputFormat);
                if iscell(counterEx)
                    break
                end
            end
        elseif isa(lb, "cell") && length(prop) == 1
            for arr = 1:length(lb)
                counterEx = falsify_single(net, lb{arr}, ub{arr}, inputSize, config.nRand, prop{1}.Hg, needReshape, inputFormat);
                if iscell(counterEx)
                    break
                end
            end
        end

        if iscell(counterEx)
            status = 0;
        else
            % Reachability
            if ~isa(lb, "cell") && length(prop) == 1
                if ~nnz(lb-ub)
                    status = 1;
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
        end

        tTime = toc(t);

        result.statusCode = status;
        if status == 0
            result.status = 'sat';
        elseif status == 1
            result.status = 'unsat';
        else
            result.status = 'unknown';
        end
        result.time = tTime;

        fprintf('  Instance %d: %s (%.1fs)\n', idx, result.status, tTime);

    catch ME
        result.statusCode = -1;
        result.status = 'error';
        result.error = ME.message;
        result.time = toc(t);

        fprintf('  Instance %d: error - %s\n', idx, ME.message);
    end

    results{i} = result;
end

totalTime = toc(totalStartTime);

%% Results summary
disp(' ');
satCount = sum(cellfun(@(r) strcmp(r.status, 'sat'), results));
unsatCount = sum(cellfun(@(r) strcmp(r.status, 'unsat'), results));
unknownCount = sum(cellfun(@(r) strcmp(r.status, 'unknown'), results));
errorCount = sum(cellfun(@(r) strcmp(r.status, 'error'), results));

disp(['ProbVer complete: ' num2str(unsatCount) ' verified, ' num2str(satCount) ' violated, ' ...
      num2str(unknownCount) ' unknown, ' num2str(errorCount) ' errors (' num2str(totalTime, '%.1f') 's)']);

% Save to CSV
fid = fopen(resultsFile, 'w');
fprintf(fid, 'index,onnx,vnnlib,status,time,error\n');
for i = 1:length(results)
    r = results{i};
    errorMsg = strrep(r.error, ',', ';');
    errorMsg = strrep(errorMsg, newline, ' ');
    fprintf(fid, '%d,%s,%s,%s,%.4f,%s\n', r.index, r.onnx, r.vnnlib, r.status, r.time, errorMsg);
end
fclose(fid);
disp(['Saved: ' resultsFile]);

end

%% ========== HELPER FUNCTIONS ==========

function IS = create_input_set(lb, ub, inputSize, needReshape)
    if ~isscalar(inputSize)
        lb = reshape(lb, inputSize);
        ub = reshape(ub, inputSize);
    end

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

    IS = ImageStar(lb, ub);

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
            if isa(net.Layers(1), 'nnet.cnn.layer.ImageInputLayer')
                xRand = dlarray(xRand, "SSCB");
            elseif isa(net.Layers(1), 'nnet.cnn.layer.FeatureInputLayer') || isa(net.Layers(1), 'nnet.onnx.layer.FeatureInputLayer')
                xRand = dlarray(xRand, "CB");
            else
                disp(['Warning: Unknown input layer type: ' class(net.Layers(1))]);
                xRand = dlarray(xRand, "SSCB");
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
        catch
            % Continue on prediction errors
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
