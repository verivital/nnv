function pairs = vnncomp_pick_instances(inst_csv, bench_root, sub, run_which)
% vnncomp_pick_instances  Resolve (onnx, vnnlib) pairs from an instances.csv.
%
%   Extracted (behavior-preserving) copy of the local `pick_instances` +
%   `ensure_decompressed` helpers in run_all_benchmarks.m, so the CSV-parsing /
%   instance-selection logic can be unit-tested without spinning up a parpool or
%   importing a network. Adopting it as the single source of truth in
%   run_all_benchmarks.m is a recommended follow-up; until then keep this file
%   equivalent to those local functions.
%
%   Reads instances.csv (3 cols: onnx_rel, vnnlib_rel, timeout_s), strips a
%   leading './' from each relative path, and returns a cell array of resolvable
%   structs {onnx_rel, vnnlib_rel, onnx, vnnlib} (decompressing .gz on demand).
%     run_which = 'first' -> the FIRST resolvable instance only.
%     run_which = 'all'   -> EVERY resolvable instance.

pairs = {};
T = readtable(inst_csv, 'Delimiter', ',', 'ReadVariableNames', false, ...
    'TextType', 'string', 'Format', '%s%s%s');
for r = 1:height(T)
    onnx_rel   = strrep(strtrim(char(T{r,1})), './', '');
    vnnlib_rel = strrep(strtrim(char(T{r,2})), './', '');
    onnx_p   = ensure_decompressed(fullfile(bench_root, sub, onnx_rel));
    vnnlib_p = ensure_decompressed(fullfile(bench_root, sub, vnnlib_rel));
    if ~isempty(onnx_p) && ~isempty(vnnlib_p)
        pairs{end+1} = struct('onnx_rel', onnx_rel, 'vnnlib_rel', vnnlib_rel, ...
            'onnx', onnx_p, 'vnnlib', vnnlib_p); %#ok<AGROW>
        if strcmp(run_which, 'first'), return; end
    end
end
end

function p = ensure_decompressed(p)
% Returns a path to a decompressed file, decompressing the .gz form if needed.
% Returns '' if neither exists.
if isfile(p), return; end
gz = [p '.gz'];
if isfile(gz)
    try
        names = gunzip(gz);   % returns the cellstr of extracted file path(s)
    catch
        p = '';
        return;
    end
    if isfile(p), return; end
    % gunzip may resolve to a different path than [p without .gz]; trust its return.
    if ~isempty(names) && isfile(names{1})
        p = names{1};
        return;
    end
end
p = '';
end
