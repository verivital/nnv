%% Master VNN-COMP 2025 runner
%
% For each benchmark folder under <bench_root>, reads its `instances.csv`
% to obtain a real (onnx, vnnlib) pair, decompresses the .gz copies on
% demand, then runs `run_vnncomp_instance` with a per-instance wall-clock
% timeout. Captures real error messages, writes a CSV log + markdown
% summary.
%
% Usage:
%
%   % default benchmark root: <repo>/vnncomp2025_benchmarks/benchmarks
%   run_all_benchmarks
%
%   % or explicit
%   run_all_benchmarks('C:\path\to\benchmarks', 120)
%
% Status codes:
%    0 = sat                    1 = unsat              2 = unknown
%   -1 = error / timeout       -2 = file missing      -3 = decompress failed

function run_all_benchmarks(varargin)

% --- Parse args ---
script_dir = fileparts(mfilename('fullpath'));
default_root = fullfile(script_dir, '..', '..', '..', '..', '..', '..', 'vnncomp2025_benchmarks', 'benchmarks');
if nargin >= 1 && ~isempty(varargin{1})
    bench_root = char(varargin{1});
else
    bench_root = char(default_root);
end
if nargin >= 2 && ~isempty(varargin{2})
    timeout_s = varargin{2};
else
    timeout_s = 120;
end
% Optional 3rd arg: cell array / string array of folder names to include
% (e.g. {'cifar100_2024','vggnet16_2022'}); empty = all folders.
if nargin >= 3 && ~isempty(varargin{3})
    only_folders = cellstr(string(varargin{3}));
else
    only_folders = {};
end
% Optional 4th arg: 'first' (one representative instance per folder, default) or
% 'all' (every resolvable instance in instances.csv -- the full competition set,
% needed for a true per-instance comparison to the VNN-COMP 2025 results).
if nargin >= 4 && ~isempty(varargin{4})
    run_which = lower(char(string(varargin{4})));
else
    run_which = 'first';
end

if ~isfolder(bench_root)
    error('Benchmark root not found: %s', bench_root);
end
fprintf('Benchmark root: %s\n', bench_root);
fprintf('Per-instance timeout: %d s\n\n', timeout_s);

addpath(script_dir);

% --- Map subfolder name -> NNV category string used by run_vnncomp_instance ---
% (the dispatcher matches via `contains`, so partial matches are OK)
subfolder_to_category = containers.Map();
subfolder_to_category('acasxu_2023')                       = 'acasxu';
subfolder_to_category('cctsdb_yolo_2023')                  = 'cctsdb_yolo';
subfolder_to_category('cersyve')                           = 'cersyve';
subfolder_to_category('cgan_2023')                         = 'cgan';
subfolder_to_category('cifar100_2024')                     = 'cifar100';
subfolder_to_category('collins_aerospace_benchmark')       = 'collins_aerospace_benchmark';
subfolder_to_category('collins_rul_cnn_2022')              = 'collins_rul';
subfolder_to_category('cora_2024')                         = 'cora';
subfolder_to_category('dist_shift_2023')                   = 'dist_shift';
subfolder_to_category('linearizenn_2024')                  = 'linearizenn';
subfolder_to_category('lsnc_relu')                         = 'lsnc_relu';
subfolder_to_category('malbeware')                         = 'malbeware';
subfolder_to_category('metaroom_2023')                     = 'metaroom';
subfolder_to_category('ml4acopf_2024')                     = 'ml4acopf';
subfolder_to_category('nn4sys')                            = 'nn4sys';
subfolder_to_category('relusplitter')                      = 'relusplitter';
subfolder_to_category('safenlp_2024')                      = 'safenlp';
subfolder_to_category('sat_relu')                          = 'sat_relu';
subfolder_to_category('soundnessbench')                    = 'soundnessbench';
subfolder_to_category('test')                              = 'test';   % no NNV dispatch; will error
subfolder_to_category('tinyimagenet_2024')                 = 'tinyimagenet';
subfolder_to_category('tllverifybench_2023')               = 'tllverifybench';
subfolder_to_category('traffic_signs_recognition_2023')    = 'traffic';
subfolder_to_category('vggnet16_2022')                     = 'vggnet';
subfolder_to_category('vit_2023')                          = 'vit';
subfolder_to_category('yolo_2023')                         = 'yolo';

% --- Discover benchmark folders ---
sub_dirs = dir(bench_root);
sub_dirs = sub_dirs([sub_dirs.isdir] & ~startsWith({sub_dirs.name}, '.'));
folders = {sub_dirs.name};
if ~isempty(only_folders)
    folders = folders(ismember(folders, only_folders));
end
n = numel(folders);

% --- Set up output ---
ts = datestr(now, 'yyyymmdd_HHMMSS');
csv_path = fullfile(script_dir, sprintf('results_%s.csv', ts));
md_path  = fullfile(script_dir, sprintf('results_%s.md', ts));
fid = fopen(csv_path, 'w');
fprintf(fid, 'subfolder,category,onnx,vnnlib,status,status_str,time_s,error_message\n');

results = cell(0, 8);   % grows per instance (all-instances mode -> count unknown up front)

% --- Use Processes pool with 1 worker for clean timeouts and real error msgs ---
existing_pool = gcp('nocreate');
if ~isempty(existing_pool), delete(existing_pool); end
pool = parpool('Processes', 1);   % 1 worker avoids the 16-vs-N config error

% --- Iterate ---
for i = 1:n
    sub = folders{i};
    if ~isKey(subfolder_to_category, sub)
        % Default the category to the FOLDER NAME -- this mirrors the competition
        % harness (which passes the benchmark folder name as the category) and
        % run_vnncomp_instance routes via contains(category, "<key>"), so e.g.
        % "cgan2026"->"cgan", "soundnessbench_2026"->"soundnessbench",
        % "relusplitter_2026"->"relusplitter" all dispatch correctly. Using '?'
        % (the old default) left every new 2026 folder unrouted -> forced error.
        cat = sub;
    else
        cat = subfolder_to_category(sub);
    end
    fprintf('\n[%2d/%d] %s -> dispatcher category "%s"\n', i, n, sub, cat);

    inst_csv = fullfile(bench_root, sub, 'instances.csv');
    if ~isfile(inst_csv)
        % VNN-COMP 2026 splits each benchmark by VNN-LIB version: <sub>/1.0/instances.csv
        % (regular track) + <sub>/2.0/ (extended/VNN-LIB-2.0 track). Prefer 1.0; fall back
        % to 2.0 so the extended-ONLY benchmarks (adaptive_cruise/isomorphic/monotonic/
        % smart_turn -- which ship 2.0/ only) are run too. The 2.0 specs dispatch through
        % load_vnnlib2 (sound-or-unknown: multi-net/nonlinear/multimodal -> unknown).
        alt = fullfile(bench_root, sub, '1.0', 'instances.csv');
        if isfile(alt)
            inst_csv = alt;
        else
            alt2 = fullfile(bench_root, sub, '2.0', 'instances.csv');
            if isfile(alt2), inst_csv = alt2; end
        end
    end
    if ~isfile(inst_csv)
        row = {sub, cat, '', '', -2, 'no_instances_csv', 0, ''};
        results(end+1,:) = row; write_row(fid, row); %#ok<AGROW>
        fprintf('  -> no instances.csv\n');
        continue;
    end
    pairs = pick_instances(inst_csv, bench_root, sub, run_which);
    if isempty(pairs)
        row = {sub, cat, '', '', -2, 'no_resolvable_instance', 0, ''};
        results(end+1,:) = row; write_row(fid, row); %#ok<AGROW>
        fprintf('  -> no resolvable instance in CSV\n');
        continue;
    end
    fprintf('  %d instance(s) [%s]\n', numel(pairs), run_which);
    % A FRESH Processes worker pays a large one-time cost on its first instance
    % (worker init + ONNX-importer custom-layer codegen): measured ~100 s even for
    % a model that solves in 9 s warm. Grant that first instance a one-time grace
    % so cold-start cannot masquerade as a solver timeout. Applies to the job's
    % first instance and to the first instance after every timeout-triggered pool
    % restart. The recorded time_s is still the true wall time.
    cold_grace_s = 120;
    pool_is_cold = true;
    for k = 1:numel(pairs)
        pr = pairs{k};
        out_file = [tempname '.txt'];
        t0 = tic;
        f = parfeval(pool, @run_vnncomp_instance, 2, cat, pr.onnx, pr.vnnlib, out_file);
        ok = wait(f, 'finished', timeout_s + cold_grace_s*pool_is_cold);
        pool_is_cold = false;
        if ok
            if isempty(f.Error)
                try
                    [status, ~] = fetchOutputs(f);
                    status_str = status_to_str(status);
                    err = '';
                catch ME
                    status = -1; status_str = 'error'; err = ME.message;
                end
            else
                status = -1; status_str = 'error'; err = f.Error.message;
            end
        else
            cancel(f);
            status = -1; status_str = 'timeout';
            err = sprintf('exceeded %d s', timeout_s);
            % cancel() CANNOT interrupt a BUSY worker (it only de-schedules queued
            % futures): the zombie keeps grinding the timed-out reach while the next
            % parfeval queues BEHIND it -- cascading false timeouts plus unbounded
            % CPU/memory growth. This is what killed the malbeware sweep runners
            % (~50 scaled_16 instances need ~370 s each vs the 120 s cap; both sweep
            % runs died with NO logs and NO rows = runner-level death). Restarting
            % the pool is the only reliable way to kill a busy MATLAB worker; the
            % ~15 s restart cost applies only on timeouts.
            try
                delete(pool);
            catch
            end
            pool = parpool('Processes', 1);
            pool_is_cold = true;   % next instance gets the one-time cold-start grace
        end
        time_s = toc(t0);
        fprintf('  [%d/%d] %s -> %s (%.1f s)%s\n', k, numel(pairs), pr.onnx_rel, status_str, time_s, ...
            ternary(~isempty(err), [': ' err(1:min(120,end))], ''));
        row = {sub, cat, pr.onnx_rel, pr.vnnlib_rel, status, status_str, time_s, err};
        results(end+1,:) = row; write_row(fid, row); %#ok<AGROW>
        if exist(out_file, 'file'), delete(out_file); end   % avoid leaving 1000s of temp .txt
    end
end

fclose(fid);
delete(pool);
fprintf('\nCSV written: %s\n', csv_path);

% --- Generate summary markdown ---
write_summary_md(md_path, results, bench_root, timeout_s);
fprintf('Summary written: %s\n', md_path);

end

%% --- Helpers ---

function pairs = pick_instances(inst_csv, bench_root, sub, run_which)
% Read instances.csv (3 cols: onnx_rel, vnnlib_rel, timeout_s) and return a cell
% array of resolvable {onnx_rel, vnnlib_rel, onnx, vnnlib} structs (decompressing
% .gz on demand). run_which='first' returns the first resolvable instance; 'all'
% returns EVERY resolvable instance (the full competition set for that benchmark).
pairs = {};
% Resolve onnx/vnnlib relative to the instances.csv's OWN directory so both the 2025
% layout (<sub>/) and the 2026 version layout (<sub>/1.0/) work (bench_root/sub then unused).
bench_dir = fileparts(inst_csv);
T = readtable(inst_csv, 'Delimiter', ',', 'ReadVariableNames', false, ...
    'TextType', 'string', 'Format', '%s%s%s');
for r = 1:height(T)
    onnx_rel   = strrep(strtrim(char(T{r,1})), './', '');
    vnnlib_rel = strrep(strtrim(char(T{r,2})), './', '');
    onnx_p   = ensure_decompressed(fullfile(bench_dir, onnx_rel));
    vnnlib_p = ensure_decompressed(fullfile(bench_dir, vnnlib_rel));
    if ~isempty(onnx_p) && ~isempty(vnnlib_p)
        pairs{end+1} = struct('onnx_rel', onnx_rel, 'vnnlib_rel', vnnlib_rel, ...
            'onnx', onnx_p, 'vnnlib', vnnlib_p); %#ok<AGROW>
        if strcmp(run_which, 'first'), return; end
    end
end
end

function p = ensure_decompressed(p)
% Returns a path to a decompressed file, decompressing the .gz form if
% needed. Returns '' if neither exists.
if isfile(p), return; end
gz = [p '.gz'];
if isfile(gz)
    try
        gunzip(gz);
    catch
        p = '';
        return;
    end
    if isfile(p), return; end
end
p = '';
end

function s = status_to_str(code)
switch code
    case 0,  s = 'sat';
    case 1,  s = 'unsat';
    case 2,  s = 'unknown';
    case -1, s = 'error';
    case -2, s = 'missing';
    case -3, s = 'decompress_failed';
    otherwise, s = sprintf('code_%d', code);
end
end

function r = ternary(c, a, b), if c, r = a; else, r = b; end, end

function write_row(fid, row)
% row = {sub, cat, onnx_rel, vnnlib_rel, status, status_str, time_s, err}
err_csv = strrep(row{8}, ',', ';');
err_csv = strrep(err_csv, sprintf('\n'), ' | ');
err_csv = strrep(err_csv, '"', '''');
fprintf(fid, '%s,%s,%s,%s,%d,%s,%.2f,"%s"\n', ...
    row{1}, row{2}, row{3}, row{4}, row{5}, row{6}, row{7}, err_csv);
end

function write_summary_md(path, results, bench_root, timeout_s)
n = size(results, 1);
fid = fopen(path, 'w');
fprintf(fid, '# VNN-COMP 2025 sweep results\n\n');
% Sanitize the bench root so committed result summaries don't leak a local
% home/Dropbox absolute path (show it repo-relative from vnncomp2025_benchmarks).
bench_root_disp = regexprep(strrep(bench_root,'\','/'), '.*/(vnncomp2025_benchmarks/.*)$', '<repo>/$1');
fprintf(fid, 'Bench root: `%s`  \n', bench_root_disp);
fprintf(fid, 'Per-instance timeout: %d s  \n', timeout_s);
fprintf(fid, 'Run at: %s  \n\n', datestr(now,'yyyy-mm-dd HH:MM:SS'));

% Tally
counts = struct('sat',0,'unsat',0,'unknown',0,'error',0,'missing',0,'timeout',0,'decompress_failed',0,'no_instances_csv',0,'no_resolvable_instance',0,'other',0);
for i = 1:n
    s = matlab.lang.makeValidName(char(results{i,6}));
    if isfield(counts, s), counts.(s) = counts.(s)+1; else, counts.other = counts.other+1; end
end
fprintf(fid, '## Tally (per recorded instance)\n\n');
fprintf(fid, '| outcome | count |\n|---|---|\n');
flds = fieldnames(counts);
for i = 1:numel(flds)
    if counts.(flds{i}) > 0
        fprintf(fid, '| %s | %d |\n', flds{i}, counts.(flds{i}));
    end
end

fprintf(fid, '\n## Per-folder result\n\n');
fprintf(fid, '| folder | category | status | time (s) | error |\n|---|---|---|---|---|\n');
for i = 1:n
    sub = results{i,1}; cat = results{i,2}; ss = results{i,6}; t = results{i,7}; err = results{i,8};
    err_md = strrep(err, '|', '\\|');
    err_md = strrep(err_md, sprintf('\n'), ' ');
    if numel(err_md) > 120, err_md = [err_md(1:120) '...']; end
    fprintf(fid, '| %s | %s | %s | %.1f | %s |\n', sub, cat, ss, t, err_md);
end

fprintf(fid, '\n## Folders grouped by outcome\n\n');
groups = struct();
for i = 1:n
    s = matlab.lang.makeValidName(char(results{i,6}));
    if ~isfield(groups, s), groups.(s) = {}; end
    groups.(s){end+1} = results{i,1};
end
flds = fieldnames(groups);
for i = 1:numel(flds)
    fprintf(fid, '- **%s** (%d): %s\n', flds{i}, numel(groups.(flds{i})), strjoin(groups.(flds{i}), ', '));
end
fclose(fid);
end
