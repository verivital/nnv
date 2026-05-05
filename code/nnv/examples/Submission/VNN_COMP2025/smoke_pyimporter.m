function smoke_pyimporter(varargin)
%SMOKE_PYIMPORTER  Run the Python ONNX-to-NNV importer + cross-validation
% across all VNN-COMP 2025 benchmark folders. For each folder:
%   1. Read instances.csv to pick first instance
%   2. Decompress the ONNX
%   3. Call `python tools/onnx2nnv_python/onnx2nnv.py` to produce .nnv.mat
%   4. Call load_nnv_from_mat to construct the NNV NN object
%   5. Cross-validate against onnxruntime on N random inputs
%   6. (Optional) try NNV reach() with approx-star to confirm it runs
%
% Output: results_pyimporter_<ts>.csv

script_dir = fileparts(mfilename('fullpath'));
default_root = fullfile(script_dir, '..', '..', '..', '..', '..', '..', 'vnncomp2025_benchmarks', 'benchmarks');
if nargin >= 1 && ~isempty(varargin{1}), bench_root = char(varargin{1}); else, bench_root = char(default_root); end
if nargin >= 2 && ~isempty(varargin{2}), reach_timeout = varargin{2}; else, reach_timeout = 15; end

repo_root = fullfile(script_dir, '..', '..', '..', '..', '..', '..');
py_importer = fullfile(repo_root, 'tools', 'onnx2nnv_python', 'onnx2nnv.py');

addpath(genpath(fullfile(script_dir, '..', '..', '..', 'engine')));
addpath(script_dir);

ts = datestr(now,'yyyymmdd_HHMMSS');
csv_path = fullfile(script_dir, sprintf('results_pyimporter_%s.csv', ts));
md_path  = fullfile(script_dir, sprintf('results_pyimporter_%s.md', ts));
fid = fopen(csv_path,'w');
fprintf(fid, 'subfolder,onnx,vnnlib,import_status,layers,xval_max_diff,reach_status,time_s,err\n');

sub_dirs = dir(bench_root);
sub_dirs = sub_dirs([sub_dirs.isdir] & ~startsWith({sub_dirs.name},'.'));
folders = {sub_dirs.name};

% Use a 1-worker Processes pool for reach() timeouts
existing = gcp('nocreate'); if ~isempty(existing), delete(existing); end
pool = parpool('Processes', 1);

results = cell(0,9);

for i = 1:numel(folders)
    sub = folders{i};
    fprintf('\n[%2d/%d] %s\n', i, numel(folders), sub);
    inst_csv = fullfile(bench_root, sub, 'instances.csv');
    if ~isfile(inst_csv)
        rec = {sub,'','','no_instances_csv',0,NaN,'-',0,''};
        results = [results; rec]; %#ok<AGROW>
        write_row(fid, rec);
        continue;
    end
    pair = pick_first_instance(inst_csv, bench_root, sub);
    if isempty(pair)
        rec = {sub,'','','no_resolvable_instance',0,NaN,'-',0,''};
        results = [results; rec]; %#ok<AGROW>
        write_row(fid, rec);
        continue;
    end

    onnx_p = pair.onnx;
    vnnlib_p = pair.vnnlib;
    fprintf('  onnx:   %s\n', onnx_p);
    fprintf('  vnnlib: %s\n', vnnlib_p);

    % --- Run Python importer ---
    mat_p = strrep(onnx_p, '.onnx', '.nnv.mat');
    cmd = sprintf('python "%s" "%s" "%s"', py_importer, onnx_p, mat_p);
    [s, o] = system(cmd);
    if s ~= 0 || ~isfile(mat_p)
        err = strtrim(o);
        if numel(err) > 200, err = [err(1:200) '...']; end
        rec = {sub, pair.onnx_rel, pair.vnnlib_rel, 'import_failed', 0, NaN, '-', 0, err};
        fprintf('  -> import_failed: %s\n', err);
        results = [results; rec]; %#ok<AGROW>
        write_row(fid, rec);
        continue;
    end

    % --- Load via load_nnv_from_mat ---
    try
        net = load_nnv_from_mat(mat_p);
        n_layers = numel(net.Layers);
        fprintf('  loaded %d layers\n', n_layers);
    catch ME
        rec = {sub, pair.onnx_rel, pair.vnnlib_rel, 'load_failed', 0, NaN, '-', 0, ME.message};
        fprintf('  -> load_failed: %s\n', ME.message);
        results = [results; rec]; %#ok<AGROW>
        write_row(fid, rec);
        continue;
    end

    % --- Cross-validate ---
    try
        xval = xvalidate(onnx_p, mat_p, 3);
        if isnan(xval)
            xval_str = 'xval_fail';
        elseif xval < 1e-3
            xval_str = 'xval_pass';
        elseif xval < 1e-1
            xval_str = 'xval_loose';
        else
            xval_str = 'xval_FAIL';
        end
        fprintf('  xval max diff: %.6e (%s)\n', xval, xval_str);
    catch ME
        xval = NaN; xval_str = 'xval_err';
        fprintf('  xval err: %s\n', ME.message);
    end

    % --- Try reach with approx-star (15s) ---
    reach_status = '-'; reach_time = 0;
    if ~isnan(xval) && xval < 1e-1
        try
            prop = load_vnnlib(vnnlib_p);
            if iscell(prop.lb), lb = prop.lb{1}; ub = prop.ub{1}; spec = prop.prop{1};
            else, lb = prop.lb; ub = prop.ub; spec = prop.prop{1}; end
            in_layer = net.Layers{1};
            if isa(in_layer, 'FeatureInputLayer')
                IS = Star(double(lb(:)), double(ub(:)));
            else
                sz = in_layer.InputSize;
                IS = ImageStar(reshape(double(lb), sz), reshape(double(ub), sz));
            end
            reachOpts = struct('reachMethod','approx-star');
            t0 = tic;
            f = parfeval(pool, @(net_in, IS_in, opts) net_in.reach(IS_in, opts), 1, net, IS, reachOpts);
            ok = wait(f, 'finished', reach_timeout);
            if ok
                if isempty(f.Error)
                    ySet = fetchOutputs(f);
                    status = verify_specification(ySet, prop.prop);
                    reach_status = status_to_str(status);
                else
                    reach_status = 'reach_err';
                end
            else
                cancel(f);
                reach_status = 'timeout';
            end
            reach_time = toc(t0);
            fprintf('  reach: %s (%.1fs)\n', reach_status, reach_time);
        catch ME
            reach_status = 'reach_setup_err'; reach_time = 0;
            fprintf('  reach setup err: %s\n', ME.message);
        end
    end

    rec = {sub, pair.onnx_rel, pair.vnnlib_rel, ...
        sprintf('import_ok/%s', xval_str), n_layers, xval, reach_status, reach_time, ''};
    results = [results; rec]; %#ok<AGROW>
    write_row(fid, rec);
end

fclose(fid);
delete(pool);
write_md(md_path, results);
fprintf('\nResults: %s\n  %s\n', csv_path, md_path);
end


function s = status_to_str(code)
switch code
    case 0,  s = 'sat';
    case 1,  s = 'unsat';
    case 2,  s = 'unknown';
    otherwise, s = sprintf('code_%d', code);
end
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
        pair = struct('onnx',onnx_p,'vnnlib',vnnlib_p, ...
                      'onnx_rel',onnx_rel,'vnnlib_rel',vnnlib_rel);
        return;
    end
end
end

function p = ensure_decompressed(p)
if isfile(p), return; end
gz = [p '.gz'];
if isfile(gz)
    try, out = gunzip(gz); catch, p = ''; return; end
    if iscell(out) && ~isempty(out) && isfile(out{1}), p = out{1}; return; end
    if isfile(p), return; end
end
p = '';
end

function write_row(fid, row)
err = strrep(row{9}, ',', ';'); err = strrep(err, sprintf('\n'),' | '); err = strrep(err,'"','''');
fprintf(fid, '%s,%s,%s,%s,%d,%g,%s,%.2f,"%s"\n', row{1},row{2},row{3},row{4},row{5},row{6},row{7},row{8},err);
end

function write_md(path, results)
n = size(results,1);
fid = fopen(path,'w');
fprintf(fid, '# VNN-COMP 2025 sweep — Python importer + cross-validation\n\n');
fprintf(fid, 'Generated: %s\n\n', datestr(now,'yyyy-mm-dd HH:MM:SS'));
fprintf(fid, '| folder | import | layers | xval_max_diff | reach | time | err |\n|---|---|---|---|---|---|---|\n');
for i=1:n
    err = results{i,9}; if numel(err)>120, err = [err(1:120) '...']; end
    fprintf(fid, '| %s | %s | %d | %.4g | %s | %.1f | %s |\n', ...
        results{i,1},results{i,4},results{i,5},results{i,6},results{i,7},results{i,8},err);
end
fclose(fid);
end
