function verify_one_instance(idx, benchmarkDir, outputCSV, nRand)
% verify_one_instance - Verify a single yolo_2023 ProbVer instance.
%
% Designed to be invoked once per instance from a bash orchestrator
% (run_probver.sh) so each instance runs in its own MATLAB address space.
% This isolates the per-instance memory footprint, so an OOM-killed run
% (e.g. inside Prob_reach for a high-coverage TinyYOLO property) takes
% down only one instance instead of the whole repeatability run.
%
% Inputs:
%   idx           - 1-based row of instances.csv to verify
%   benchmarkDir  - absolute path to yolo_2023 directory
%   outputCSV     - results CSV; one row appended on success/error
%   nRand         - number of random falsification samples
%
% On the OOM path MATLAB itself receives SIGKILL; the bash orchestrator
% records that as `oom` in the CSV. On the happy path this function
% appends one row and exits 0.

% Bootstrap NNV path + GPU forward-compat (Blackwell / RTX 5090).
scriptDir = fileparts(mfilename('fullpath'));
nnvRoot = fullfile(scriptDir, '..', '..', '..');
if exist(fullfile(nnvRoot, 'startup_nnv.m'), 'file')
    addpath(genpath(nnvRoot));
end
try
    parallel.gpu.enableCUDAForwardCompatibility(true);
catch
end

% Resolve the instance.
instancesFile = fullfile(benchmarkDir, 'instances.csv');
fid = fopen(instancesFile, 'r');
lines = {};
while ~feof(fid)
    line = fgetl(fid);
    if ischar(line) && ~isempty(strtrim(line))
        lines{end+1} = line; %#ok<AGROW>
    end
end
fclose(fid);

if idx < 1 || idx > numel(lines)
    error('verify_one_instance:badIndex', 'idx %d out of range 1..%d', idx, numel(lines));
end
parts = strsplit(lines{idx}, ',');
onnxRel = strtrim(parts{1});
vnnlibRel = strtrim(parts{2});

% Load network (per-process, fresh state).
onnxPath = fullfile(benchmarkDir, 'onnx', 'TinyYOLO.onnx');
net = importNetworkFromONNX(onnxPath, 'InputDataFormats', 'BCSS');
inputSize = net.Layers(1, 1).InputSize;
inputFormat = 'default';
needReshape = 2;

reachOptions = struct;
reachOptions.train_epochs       = 200;
reachOptions.train_lr           = 0.0001;
reachOptions.coverage           = 0.999;
reachOptions.confidence         = 0.999;
reachOptions.train_mode         = 'Linear';
reachOptions.surrogate_dim      = [-1, -1];
reachOptions.threshold_normal   = 1e-5;
reachOptions.reachMethod        = 'cp-star';

vnnlibPath = fullfile(benchmarkDir, vnnlibRel);
property = load_vnnlib(vnnlibPath);
lb = property.lb; ub = property.ub; prop = property.prop;

statusStr = 'unknown';
errMsg = '';
tStart = posixtime(datetime('now'));

try
    counterEx = nan;
    if ~isa(lb, 'cell') && length(prop) == 1
        counterEx = local_falsify(net, lb, ub, inputSize, nRand, prop{1}.Hg, needReshape, inputFormat);
    elseif isa(lb, 'cell') && length(lb) == length(prop)
        for spc = 1:length(lb)
            counterEx = local_falsify(net, lb{spc}, ub{spc}, inputSize, nRand, prop{spc}.Hg, needReshape, inputFormat);
            if iscell(counterEx); break; end
        end
    elseif isa(lb, 'cell') && length(prop) == 1
        for arr = 1:length(lb)
            counterEx = local_falsify(net, lb{arr}, ub{arr}, inputSize, nRand, prop{1}.Hg, needReshape, inputFormat);
            if iscell(counterEx); break; end
        end
    end

    if iscell(counterEx)
        statusStr = 'sat';
    else
        if ~isa(lb, 'cell') && length(prop) == 1
            if ~nnz(lb - ub)
                statusStr = 'unsat';
            else
                IS = local_input_set(lb, ub, inputSize, needReshape);
                ySet = Prob_reach(net, IS, reachOptions);
                statusStr = local_status_to_str(verify_specification(ySet, prop));
            end
        elseif isa(lb, 'cell') && length(lb) == length(prop)
            local_status = 2 * ones(length(lb), 1);
            for spc = 1:length(lb)
                lb_spc = lb{spc}; ub_spc = ub{spc};
                if ~nnz(lb_spc - ub_spc)
                    local_status(spc) = 1;
                else
                    IS = local_input_set(lb_spc, ub_spc, inputSize, needReshape);
                    ySet = Prob_reach(net, IS, reachOptions);
                    if isempty(ySet.C)
                        dd = ySet.V; DD = ySet.V;
                        ySet = Star(dd, DD);
                    end
                    local_status(spc) = verify_specification(ySet, prop(spc));
                end
            end
            statusStr = local_status_to_str(all(local_status == 1));
        elseif isa(lb, 'cell') && length(prop) == 1
            local_status = 2 * ones(length(lb), 1);
            for spc = 1:length(lb)
                lb_spc = lb{spc}; ub_spc = ub{spc};
                if ~nnz(lb_spc - ub_spc)
                    local_status(spc) = 1;
                else
                    IS = local_input_set(lb_spc, ub_spc, inputSize, needReshape);
                    ySet = Prob_reach(net, IS, reachOptions);
                    local_status(spc) = verify_specification(ySet, prop(spc));
                end
            end
            statusStr = local_status_to_str(all(local_status == 1));
        end
    end
catch ME
    statusStr = 'error';
    errMsg = strrep(ME.message, ',', ';');
    errMsg = strrep(errMsg, newline, ' ');
end

elapsed = posixtime(datetime('now')) - tStart;

fid = fopen(outputCSV, 'a');
fprintf(fid, '%d,%s,%s,%s,%.4f,%s\n', idx, onnxRel, vnnlibRel, statusStr, elapsed, errMsg);
fclose(fid);

fprintf('verify_one_instance: idx=%d status=%s time=%.4fs\n', idx, statusStr, elapsed);
end

% =====================================================================
% Local helpers — duplicated from run_probver.m so this file is
% self-contained.
% =====================================================================

function s = local_status_to_str(c)
    if islogical(c)
        if c, s = 'unsat'; else, s = 'unknown'; end
        return;
    end
    switch c
        case 0, s = 'sat';
        case 1, s = 'unsat';
        otherwise, s = 'unknown';
    end
end

function IS = local_input_set(lb, ub, inputSize, needReshape)
    if ~isscalar(inputSize)
        lb = reshape(lb, inputSize);
        ub = reshape(ub, inputSize);
    end
    if needReshape == 1
        lb = permute(lb, [2 1 3]);
        ub = permute(ub, [2 1 3]);
    elseif needReshape == 2
        newSize = [inputSize(2) inputSize(1) inputSize(3:end)];
        lb = reshape(lb, newSize); lb = permute(lb, [2 1 3 4]);
        ub = reshape(ub, newSize); ub = permute(ub, [2 1 3 4]);
    end
    IS = ImageStar(lb, ub);
    try
        xxx = find((lb - ub));
        xxx = reshape(xxx, [], 1);
        if numel(lb) ~= length(xxx)
            IS.C = IS.C(:, xxx);
            IS.pred_lb = IS.pred_lb(xxx);
            IS.pred_ub = IS.pred_ub(xxx);
            xxx = xxx + 1;
            IS.V = IS.V(:, :, :, [1; xxx]);
            IS.numPred = length(xxx);
        end
    end
end

function xRand = local_random_examples(net, lb, ub, nR, inputSize, needReshape, inputFormat)
    xB = Box(lb, ub);
    xRand = xB.sample(nR - 2);
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
        if strcmp(inputFormat, 'default')
            if isa(net.Layers(1, 1), 'nnet.cnn.layer.ImageInputLayer')
                xRand = dlarray(xRand, 'SSCB');
            elseif isa(net.Layers(1, 1), 'nnet.cnn.layer.FeatureInputLayer') || ...
                   isa(net.Layers(1, 1), 'nnet.onnx.layer.FeatureInputLayer')
                xRand = dlarray(xRand, 'CB');
            else
                error('verify_one_instance:unknownInputLayer', 'Unknown input format');
            end
        else
            if contains(inputFormat, 'U')
                xRand = dlarray(xRand, inputFormat + 'U');
            else
                xRand = dlarray(xRand, inputFormat);
            end
        end
    end
end

function counterEx = local_falsify(net, lb, ub, inputSize, nRand, Hs, needReshape, inputFormat)
    counterEx = nan;
    xRand = local_random_examples(net, lb, ub, nRand, inputSize, needReshape, inputFormat);
    s = size(xRand);
    n = length(s);
    for i = 1:s(n)
        x = local_get_example(xRand, i);
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
        end
    end
end

function x = local_get_example(xRand, i)
    s = size(xRand);
    n = length(s);
    if n == 4
        x = xRand(:, :, :, i);
    elseif n == 3
        x = xRand(:, :, i);
    elseif n == 2
        x = xRand(:, i);
        xsize = size(x);
        if xsize(1) ~= 1 && ~isa(x, 'dlarray')
            x = x';
        end
    else
        error('verify_one_instance:badInput', 'InputSize = %s', mat2str(s));
    end
end
