%% Smoke test: MATLAB built-in CROWN across all 26 VNN-COMP 2025 categories
%
% For each benchmark folder, try a sequence of strategies to obtain a clean
% dlnetwork (featureInputLayer / imageInputLayer + built-in body layers),
% then call verifyNetworkRobustness(net, lb, ub, target). Records which
% categories MATLAB's built-in verifier can handle.
%
% Records to results_matlab_crown_<ts>.csv with columns:
%   subfolder,strategy,n_layers,n_custom,verify_status,time_s,err
%
% Strategies tried per category (first that gives an all-built-in net wins):
%   S1. importNetworkFromONNX with various InputDataFormats hints
%   S2. importNetworkFromONNX + rebuild as featureInputLayer body-only
%
% Status:
%   verified  = formal proof of robustness (MATLAB returns 'verified')
%   violated  = counterexample / proven not-robust
%   unproven  = bound too loose to decide
%   import_failed
%   custom_layer (importer left an opaque non-built-in layer)
%   verify_failed
%   no_files

function smoke_matlab_crown(varargin)

script_dir = fileparts(mfilename('fullpath'));
default_root = fullfile(script_dir, '..', '..', '..', '..', '..', '..', 'vnncomp2025_benchmarks', 'benchmarks');
if nargin >= 1 && ~isempty(varargin{1}), bench_root = char(varargin{1}); else, bench_root = char(default_root); end
if nargin >= 2 && ~isempty(varargin{2}), per_inst_timeout = varargin{2}; else, per_inst_timeout = 60; end

addpath('c:/Users/taylo/Dropbox/Research/talks/vu-isis-2025-11-15/vsc_matlab_mcp/nnv/code/nnv');

ts = datestr(now,'yyyymmdd_HHMMSS');
csv_path = fullfile(script_dir, sprintf('results_matlab_crown_%s.csv', ts));
md_path  = fullfile(script_dir, sprintf('results_matlab_crown_%s.md', ts));
fid = fopen(csv_path,'w');
fprintf(fid, 'subfolder,strategy,n_layers,n_custom,status,time_s,err\n');

sub_dirs = dir(bench_root);
sub_dirs = sub_dirs([sub_dirs.isdir] & ~startsWith({sub_dirs.name},'.'));
folders = {sub_dirs.name};

results = cell(0,7);

for i = 1:numel(folders)
    sub = folders{i};
    fprintf('\n[%2d/%d] %s\n', i, numel(folders), sub);
    inst_csv = fullfile(bench_root, sub, 'instances.csv');
    if ~isfile(inst_csv), append_row(fid, results, {sub,'-',0,0,'no_instances_csv',0,''}); continue; end
    pair = pick_first_instance(inst_csv, bench_root, sub);
    if isempty(pair), append_row(fid, results, {sub,'-',0,0,'no_files',0,''}); continue; end

    % --- Strategy 1: try importNetworkFromONNX with format hints ---
    net = []; n_custom = inf; chosen_fmt = '';
    fmts = {{'BCSS','BC'}, {'BC','BC'}, {'BCT','BC'}, {'BSSC','BC'}, {'TBC','BC'}, {'BTC','BC'}};
    warning('off','all');
    for fi = 1:numel(fmts)
        fmt = fmts{fi};
        try
            try_net = importNetworkFromONNX(pair.onnx, 'InputDataFormats', fmt{1}, 'OutputDataFormats', fmt{2});
            n_c = count_custom(try_net);
            if n_c < n_custom, net = try_net; n_custom = n_c; chosen_fmt = strjoin(fmt,'/'); end
            if n_c == 0, break; end
        catch
            % skip
        end
    end
    warning('on','all');

    if isempty(net)
        append_row(fid, results, {sub,'import',0,inf,'import_failed',0,''});
        continue;
    end
    n_layers = numel(net.Layers);
    if n_custom > 0
        % Strategy 2: try to rebuild as plain featureInputLayer + extracted FC/Relu
        rebuilt = try_rebuild_feature(net);
        if ~isempty(rebuilt)
            net = rebuilt; n_custom = 0;
        else
            append_row(fid, results, {sub,sprintf('fmt=%s',chosen_fmt),n_layers,n_custom,'custom_layer',0,''});
            continue;
        end
    end

    % --- Run CROWN ---
    prop = load_vnnlib(pair.vnnlib);
    if iscell(prop.lb), lb_v = prop.lb{1}; ub_v = prop.ub{1};
    else, lb_v = prop.lb; ub_v = prop.ub; end
    lb_v = single(lb_v); ub_v = single(ub_v);

    try
        [lb_dl, ub_dl, target] = build_inputs_and_target(net, lb_v, ub_v);
    catch ME
        append_row(fid, results, {sub,sprintf('fmt=%s',chosen_fmt),numel(net.Layers),0,'input_format_failed',0,ME.message});
        continue;
    end

    t0 = tic;
    try
        result = verifyNetworkRobustness(net, lb_dl, ub_dl, target);
        status = char(string(result));
        elapsed = toc(t0);
        fprintf('  -> %s (%.1f s, %d layers, %d custom)\n', status, elapsed, numel(net.Layers), n_custom);
        append_row(fid, results, {sub,sprintf('fmt=%s',chosen_fmt),numel(net.Layers),n_custom,status,elapsed,''});
    catch ME
        elapsed = toc(t0);
        append_row(fid, results, {sub,sprintf('fmt=%s',chosen_fmt),numel(net.Layers),n_custom,'verify_failed',elapsed,ME.message});
        fprintf('  -> verify_failed: %s\n', ME.message);
    end
end

fclose(fid);
write_md(md_path, results);
fprintf('\nResults: %s\n%s\n', csv_path, md_path);

end


function n = count_custom(net)
    n = 0;
    for j = 1:numel(net.Layers)
        cls = class(net.Layers(j));
        if ~startsWith(cls,'nnet.cnn.layer.') && ~startsWith(cls,'nnet.onnx.layer.') ...
                && ~startsWith(cls,'nnet.keras.layer.')
            n = n + 1;
        end
    end
end

function net2 = try_rebuild_feature(net)
% Try to build a featureInputLayer + FC/ReLU body, skipping reshape/flatten/permute custom layers
    body = {};
    in_dim = [];
    for j = 1:numel(net.Layers)
        L = net.Layers(j);
        if isa(L,'nnet.cnn.layer.FeatureInputLayer')
            in_dim = L.InputSize; if iscolumn(in_dim), in_dim = in_dim(1); end
        elseif isa(L,'nnet.cnn.layer.ImageInputLayer')
            in_dim = prod(L.InputSize);
        elseif isa(L,'nnet.cnn.layer.SequenceInputLayer')
            in_dim = L.InputSize; if iscolumn(in_dim), in_dim = in_dim(1); end
        elseif isa(L,'nnet.cnn.layer.FullyConnectedLayer')
            body{end+1} = fullyConnectedLayer(L.OutputSize,'Name',L.Name,'Weights',L.Weights,'Bias',L.Bias); %#ok<AGROW>
        elseif isa(L,'nnet.cnn.layer.ReLULayer')
            body{end+1} = reluLayer('Name',L.Name); %#ok<AGROW>
        elseif isa(L,'nnet.cnn.layer.LeakyReLULayer')
            body{end+1} = leakyReluLayer(L.Scale,'Name',L.Name); %#ok<AGROW>
        elseif isa(L,'nnet.cnn.layer.SigmoidLayer')
            body{end+1} = sigmoidLayer('Name',L.Name); %#ok<AGROW>
        elseif isa(L,'nnet.cnn.layer.TanhLayer')
            body{end+1} = tanhLayer('Name',L.Name); %#ok<AGROW>
        else
            % unknown layer breaks the rebuild
            net2 = []; return;
        end
    end
    if isempty(body) || isempty(in_dim), net2 = []; return; end
    layers = [featureInputLayer(in_dim,'Name','input'); body{:}];
    try
        net2 = dlnetwork(layers);
    catch
        net2 = [];
    end
end

function [lb_dl, ub_dl, target] = build_inputs_and_target(net, lb_v, ub_v)
    in_layer = net.Layers(1);
    if isa(in_layer,'nnet.cnn.layer.FeatureInputLayer')
        lb_dl = dlarray(lb_v(:),'CB');
        ub_dl = dlarray(ub_v(:),'CB');
        center = (lb_v(:)+ub_v(:))/2;
        yC = predict(net, dlarray(center,'CB'));
    elseif isa(in_layer,'nnet.cnn.layer.ImageInputLayer')
        sz = in_layer.InputSize;
        lb_dl = dlarray(reshape(lb_v,[sz 1]),'SSCB');
        ub_dl = dlarray(reshape(ub_v,[sz 1]),'SSCB');
        center = (lb_v+ub_v)/2;
        yC = predict(net, dlarray(reshape(center,[sz 1]),'SSCB'));
    else
        error('unsupported input layer type: %s', class(in_layer));
    end
    yC = extractdata(squeeze(yC));
    [~, target] = max(yC(:));
end

function append_row(fid, results, row)
    err = strrep(row{7}, ',', ';'); err = strrep(err, sprintf('\n'), ' | '); err = strrep(err,'"','''');
    fprintf(fid,'%s,%s,%d,%d,%s,%.2f,"%s"\n', row{1},row{2},row{3},row{4},row{5},row{6},err);
end

function write_md(path, results)
    n = size(results,1);
    fid = fopen(path,'w');
    fprintf(fid,'# MATLAB built-in CROWN smoke test (VNN-COMP 2025)\n\n');
    fprintf(fid,'| folder | strategy | layers | custom | status | time(s) | err |\n|---|---|---|---|---|---|---|\n');
    for i=1:n
        err = results{i,7}; if numel(err)>120, err = [err(1:120) '...']; end
        fprintf(fid,'| %s | %s | %d | %d | %s | %.1f | %s |\n', results{i,1},results{i,2},results{i,3},results{i,4},results{i,5},results{i,6},err);
    end
    fclose(fid);
end

function pair = pick_first_instance(inst_csv, bench_root, sub)
pair = [];
T = readtable(inst_csv,'Delimiter',',','ReadVariableNames',false,'TextType','string','Format','%s%s%s');
for r = 1:height(T)
    onnx_rel = strrep(strtrim(char(T{r,1})),'./','');
    vnnlib_rel = strrep(strtrim(char(T{r,2})),'./','');
    onnx_p = fullfile(bench_root,sub,onnx_rel);
    vnnlib_p = fullfile(bench_root,sub,vnnlib_rel);
    onnx_p = ensure_decompressed(onnx_p);
    vnnlib_p = ensure_decompressed(vnnlib_p);
    if ~isempty(onnx_p) && ~isempty(vnnlib_p)
        pair = struct('onnx',onnx_p,'vnnlib',vnnlib_p);
        return;
    end
end
end

function p = ensure_decompressed(p)
if isfile(p), return; end
gz = [p '.gz'];
if isfile(gz)
    try, gunzip(gz); catch, p = ''; return; end
    if isfile(p), return; end
end
p = '';
end
