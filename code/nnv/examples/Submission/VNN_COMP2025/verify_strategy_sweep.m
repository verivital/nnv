function verify_strategy_sweep(varargin)
%VERIFY_STRATEGY_SWEEP  Run each xval-passing benchmark through a sequence of
% reach methods and time budgets, measuring real verification verdicts.
%
% Strategy phases (in order):
%   1) approx-star, 60s
%   2) approx-star, 300s        (only if (1) was unknown/timeout)
%   3) relax-star R=0.25, 60s   (only if (1)+(2) were unknown)
%   4) relax-star R=0.5,  60s
%   5) relax-star R=0.75, 60s
%   6) exact-star, 60s          (skipped if NumStateNeurons too high — heuristic)
%
% Earliest "unsat" or "sat" verdict wins. Otherwise final verdict is "unknown".
%
% Usage:
%   verify_strategy_sweep                                % all benchmarks
%   verify_strategy_sweep({'acasxu_2023','cifar100_2024'}) % subset
%
% Output: results_strategy_<ts>.csv with one row per (benchmark, phase).

script_dir = fileparts(mfilename('fullpath'));
default_root = fullfile(script_dir, '..', '..', '..', '..', '..', '..', 'vnncomp2025_benchmarks', 'benchmarks');
bench_root = char(default_root);

if nargin >= 1 && ~isempty(varargin{1})
    benchmarks = varargin{1};
    if ischar(benchmarks) || isstring(benchmarks), benchmarks = cellstr(benchmarks); end
else
    sub_dirs = dir(bench_root);
    sub_dirs = sub_dirs([sub_dirs.isdir] & ~startsWith({sub_dirs.name},'.'));
    benchmarks = {sub_dirs.name};
end

addpath(genpath(fullfile(script_dir, '..', '..', '..', 'engine')));
addpath(script_dir);

ts = datestr(now,'yyyymmdd_HHMMSS');
csv_path = fullfile(script_dir, sprintf('results_strategy_%s.csv', ts));
fid = fopen(csv_path,'w');
fprintf(fid, 'benchmark,phase,method,relaxFactor,budget_s,verdict,time_s,err\n');

% Phase definitions
phases = {
    1, 'approx-star',  0,    60;
    2, 'approx-star',  0,    300;
    3, 'relax-star-area', 25, 60;
    4, 'relax-star-area', 50, 60;
    5, 'relax-star-area', 75, 60;
    6, 'exact-star',   0,    60;
};

% Use a process pool so we can hard-cancel timeouts
existing = gcp('nocreate'); if ~isempty(existing), delete(existing); end
pool = parpool('Processes', 1);

for i = 1:numel(benchmarks)
    sub = benchmarks{i};
    fprintf('\n[%d/%d] %s\n', i, numel(benchmarks), sub);
    inst = fullfile(bench_root, sub, 'instances.csv');
    if ~isfile(inst), continue; end
    pair = pick_first_instance(inst, bench_root, sub);
    if isempty(pair), continue; end
    onnx_p = pair.onnx; vnnlib_p = pair.vnnlib;

    mat_p = strrep(onnx_p, '.onnx', '.nnv.mat');
    if ~isfile(mat_p)
        % Re-import once per benchmark
        py_importer = fullfile(script_dir, '..', '..', '..', '..', '..', '..', 'tools', 'onnx2nnv_python', 'onnx2nnv.py');
        cmd = sprintf('python "%s" "%s" "%s"', py_importer, onnx_p, mat_p);
        [s, ~] = system(cmd);
        if s ~= 0, fprintf('  importer fail; skip\n'); continue; end
    end
    try
        net = load_nnv_from_mat(mat_p);
    catch ME
        fprintf('  load fail: %s\n', ME.message);
        continue;
    end
    try
        prop = load_vnnlib(vnnlib_p);
        if iscell(prop.lb), lb = prop.lb{1}; ub = prop.ub{1};
        else, lb = prop.lb; ub = prop.ub; end
        in_layer = net.Layers{1};
        if isa(in_layer, 'FeatureInputLayer')
            IS = Star(double(lb(:)), double(ub(:)));
        else
            sz = in_layer.InputSize;
            IS = ImageStar(reshape(double(lb), sz), reshape(double(ub), sz));
        end
    catch ME
        fprintf('  spec/IS setup fail: %s\n', ME.message);
        continue;
    end

    final_verdict = 'unknown';
    final_method = '';
    final_time = 0;

    for p = 1:size(phases,1)
        phase_idx = phases{p,1};
        method    = phases{p,2};
        rf        = phases{p,3};
        budget    = phases{p,4};

        % Skip exact-star for large nets
        if strcmp(method, 'exact-star') && numel(net.Layers) > 30
            fprintf('  phase %d %s skipped (large net)\n', phase_idx, method);
            fprintf(fid, '%s,%d,%s,%d,%d,skipped,0,too_large\n', sub, phase_idx, method, rf, budget);
            continue;
        end

        reachOpts = struct('reachMethod', method, 'relaxFactor', rf);
        t0 = tic;
        f = parfeval(pool, @(n,IS,o) n.reach(IS, o), 1, net, IS, reachOpts);
        ok = wait(f, 'finished', budget);
        elapsed = toc(t0);
        verdict = 'timeout'; err = '';
        if ok
            if isempty(f.Error)
                try
                    ySet = fetchOutputs(f);
                    status = verify_specification(ySet, prop.prop);
                    verdict = status_to_str(status);
                catch ME
                    verdict = 'verify_err';
                    err = ME.message;
                end
            else
                verdict = 'reach_err';
                err = f.Error.message;
            end
        else
            cancel(f);
        end

        % escape err for CSV
        err = strrep(err, ',', ';');
        if numel(err) > 200, err = err(1:200); end
        fprintf(fid, '%s,%d,%s,%d,%d,%s,%.2f,%s\n', sub, phase_idx, method, rf, budget, verdict, elapsed, err);
        fprintf('  phase %d %-18s rf=%d budget=%ds -> %s (%.1fs)\n', phase_idx, method, rf, budget, verdict, elapsed);

        if any(strcmp(verdict, {'unsat','sat'}))
            final_verdict = verdict;
            final_method  = method;
            final_time    = elapsed;
            break;
        end
    end

    fprintf('  ===> FINAL: %s (%s, %.1fs)\n', final_verdict, final_method, final_time);
end

fclose(fid);
delete(pool);
fprintf('\nResults: %s\n', csv_path);
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
    onnx_p = fullfile(bench_root, sub, onnx_rel);
    vnnlib_p = fullfile(bench_root, sub, vnnlib_rel);
    if ~isfile(onnx_p)
        gz = [onnx_p '.gz'];
        if isfile(gz), gunzip(gz); end
    end
    if ~isfile(vnnlib_p)
        gz = [vnnlib_p '.gz'];
        if isfile(gz), gunzip(gz); end
    end
    if isfile(onnx_p) && isfile(vnnlib_p)
        pair.onnx = onnx_p;
        pair.vnnlib = vnnlib_p;
        pair.onnx_rel = onnx_rel;
        pair.vnnlib_rel = vnnlib_rel;
        return;
    end
end
end
