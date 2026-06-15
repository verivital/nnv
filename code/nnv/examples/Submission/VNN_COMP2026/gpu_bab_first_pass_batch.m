function summary = gpu_bab_first_pass_batch(bench_root, subfolder, run_which, maxInst, goldCsv)
% GPU_BAB_FIRST_PASS_BATCH  Offline, READ-ONLY GPU-BaB first-pass over a benchmark's instances.
%
%   summary = GPU_BAB_FIRST_PASS_BATCH(bench_root, subfolder, run_which, maxInst, goldCsv)
%   resolves the (onnx, vnnlib) pairs with the SAME loader as the sweep
%   (vnncomp_pick_instances), runs the sound gpu_bab_first_pass on each, and (if goldCsv is a
%   sweep results_*.csv) classifies each instance vs gold:
%     AGREE     -- gpu robust&gold unsat, or gpu unsafe&gold sat
%     WIN       -- gpu robust/unsafe where gold was unknown/timeout/error (GPU-BaB verified
%                  what Star could not)
%     RED-FLAG  -- gpu robust but gold sat (a false-robust!), or gpu unsafe but gold unsat
%                  -- MUST be investigated before trusting GPU-BaB on this benchmark
%   It changes NO production result. RED-FLAG=0 is the soundness gate for this benchmark.
%
%   run_which: 'first' | 'all' (default 'all'); maxInst caps the count (default Inf).

    if nargin < 3 || isempty(run_which), run_which = 'all'; end
    if nargin < 4 || isempty(maxInst), maxInst = Inf; end
    if nargin < 5, goldCsv = ''; end
    cat = subfolder;   % category = folder name (matches run_all_benchmarks default routing)

    inst_csv = fullfile(bench_root, subfolder, 'instances.csv');
    if ~isfile(inst_csv)
        a1 = fullfile(bench_root, subfolder, '1.0', 'instances.csv'); if isfile(a1), inst_csv = a1; end
    end
    assert(isfile(inst_csv), 'no instances.csv for %s', subfolder);
    pairs = vnncomp_pick_instances(inst_csv, bench_root, subfolder, run_which);
    n = min(numel(pairs), maxInst);

    gold = containers.Map();
    if ~isempty(goldCsv) && isfile(goldCsv)
        G = readtable(goldCsv, 'TextType', 'string', 'Delimiter', ',');
        for r = 1:height(G)
            gold(char(G.onnx(r)) + "|" + char(G.vnnlib(r))) = lower(char(G.status_str(r)));
        end
    end

    outcsv = sprintf('gpubab_firstpass_%s.csv', subfolder);
    fid = fopen(outcsv, 'w');
    fprintf(fid, 'onnx,vnnlib,gpu_verdict,reason,guardErr,nodes,time_s,gold,class\n');
    cnt = struct('robust',0,'unsafe',0,'unknown',0,'skip',0,'error',0);
    cls = struct('AGREE',0,'WIN',0,'REDFLAG',0,'NA',0);
    opts = struct('precision','double','maxNodes',200,'nSample',20);

    for k = 1:n
        pr = pairs{k}; t0 = tic;
        try
            res = gpu_bab_first_pass(cat, pr.onnx, pr.vnnlib, opts);
        catch ME
            res = struct('verdict','error','reason',ME.message,'guardErr',NaN,'nodes',0);
        end
        ts = toc(t0);
        v = res.verdict; if isfield(cnt, v), cnt.(v) = cnt.(v) + 1; end

        gkey = pr.onnx_rel + "|" + pr.vnnlib_rel;
        gv = ''; if isKey(gold, char(gkey)), gv = gold(char(gkey)); end
        c = i_classify(v, gv);
        if isfield(cls, c), cls.(c) = cls.(c) + 1; end

        fprintf(fid, '%s,%s,%s,%s,%.2e,%d,%.1f,%s,%s\n', pr.onnx_rel, pr.vnnlib_rel, v, ...
            strrep(res.reason, ',', ';'), res.guardErr, res.nodes, ts, gv, c);
        fprintf('[%d/%d] %s -> gpu=%s gold=%s [%s] (%.1fs) %s\n', k, n, pr.onnx_rel, v, gv, c, ts, res.reason);
    end
    fclose(fid);

    summary = struct('subfolder',subfolder,'n',n,'counts',cnt,'classes',cls,'csv',outcsv);
    fprintf('\n=== GPU-BaB first-pass %s (%d instances) ===\n', subfolder, n);
    fprintf('  verdicts: robust=%d unsafe=%d unknown=%d skip=%d error=%d\n', ...
        cnt.robust, cnt.unsafe, cnt.unknown, cnt.skip, cnt.error);
    if ~isempty(goldCsv)
        fprintf('  vs gold: AGREE=%d  WIN(beat Star)=%d  RED-FLAG=%d  (RED-FLAG MUST be 0)\n', ...
            cls.AGREE, cls.WIN, cls.REDFLAG, cls.NA);
    end
end

function c = i_classify(gpu, gold)
% Classify a GPU-BaB verdict against the gold (Star) status string.
    c = 'NA';
    if isempty(gold), return; end
    holds = any(strcmp(gold, {'unsat','holds','robust','verified'}));
    viol  = any(strcmp(gold, {'sat','violated','falsified','unsafe'}));
    unc   = any(strcmp(gold, {'unknown','timeout','error'}));
    if strcmp(gpu, 'robust')
        if holds, c = 'AGREE'; elseif viol, c = 'REDFLAG'; elseif unc, c = 'WIN'; end
    elseif strcmp(gpu, 'unsafe')
        if viol, c = 'AGREE'; elseif holds, c = 'REDFLAG'; elseif unc, c = 'WIN'; end
    else
        c = 'NA';   % gpu unknown/skip -> no claim
    end
end
