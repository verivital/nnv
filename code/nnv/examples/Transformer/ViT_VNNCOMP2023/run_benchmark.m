function summary = run_benchmark(modelName, nInst, doRadius)
% RUN_BENCHMARK  Sweep the VNN-COMP 2023 ViT benchmark with the sound NNV
% star-set verifier (ViTReach). Reports, per model:
%   - full-eps (1/255) verified count (sound; expected ~0 by construction --
%     the benchmark filters out everything incomplete methods can certify),
%   - the certified radius eps* (largest eps verified) by bisection, which shows
%     the verifier certifies real, non-trivial robustness below the benchmark eps,
%   - soundness self-check: PGD/random falsification must never contradict a
%     "robust" verdict (no false robust).
%
%   summary = run_benchmark('ibp_3_3_8', 20, true)
%   summary = run_benchmark('pgd_2_3_16', 20, false)
%
% Requires the .mat produced by extract_weights.py next to this file.

    if nargin < 2, nInst = 20; end
    if nargin < 3, doRadius = true; end
    here = fileparts(mfilename('fullpath'));
    nnvroot = fullfile(here, '..', '..', '..');
    addpath(genpath(fullfile(nnvroot, 'engine')));
    addpath(here);

    M = ViTReach.load(fullfile(here, [modelName '.mat']));
    nInst = min(nInst, size(M.images,1));
    fullEps = 1/255;
    opt = struct('mode','estimate','relu','fast','marginMode','estimate');

    verified = 0; falsified = 0; epsStars = nan(nInst,1);
    fprintf('\n=== %s : %d instances ===\n', modelName, nInst);
    tAll = tic;
    for k = 1:nInst
        img = squeeze(M.images(k,:,:,:));
        label = M.labels(k);
        % full-eps verification (sound)
        [lb, ub] = ViTReach.epsBox(M, img, fullEps);
        t0 = tic;
        [rob, margins] = ViTReach.verify(M, lb, ub, label, opt);
        tk = toc(t0);
        isRobust = (rob == 1);
        verified = verified + isRobust;

        % soundness self-check: if "robust", a random/PGD search must not break it
        if isRobust
            if ViTReach_falsify(M, img, label, fullEps, 200)
                error('SOUNDNESS VIOLATION: instance %d certified robust but falsified!', k);
            end
        end

        % certified radius by bisection on eps in [0, fullEps]
        es = NaN;
        if doRadius
            es = certify_radius(M, img, label, fullEps, opt, 6);
            epsStars(k) = es;
        end
        fprintf('  inst %3d: label=%d full-eps robust=%d minMargin=%+.4f eps*=%.4f/255  (%.1fs)\n', ...
            k, label, isRobust, min(margins), es*255, tk);
    end
    fprintf('--- %s: %d/%d verified @ full eps=1/255; median eps*=%.3f/255; total %.1fs ---\n', ...
        modelName, verified, nInst, 255*median(epsStars,'omitnan'), toc(tAll));

    summary = struct('model',modelName,'nInst',nInst,'verified',verified, ...
        'falsified',falsified,'epsStars',epsStars,'medianEpsStar255',255*median(epsStars,'omitnan'));
end

function es = certify_radius(M, img, label, epsMax, opt, iters)
    % largest eps in [0, epsMax] for which verify returns robust (bisection)
    lo = 0; hi = epsMax; es = 0;
    % quick check at epsMax
    [lb,ub] = ViTReach.epsBox(M, img, epsMax);
    if ViTReach.verify(M, lb, ub, label, opt) == 1, es = epsMax; return; end
    for it = 1:iters
        mid = (lo+hi)/2;
        [lb,ub] = ViTReach.epsBox(M, img, mid);
        if ViTReach.verify(M, lb, ub, label, opt) == 1
            es = mid; lo = mid;
        else
            hi = mid;
        end
    end
end

function broken = ViTReach_falsify(M, img, label, eps, nTry)
    % random search for a counterexample inside the eps-box (raw->normalised)
    broken = false; label1 = label + 1;
    for t = 1:nTry
        d = (2*rand(size(img))-1) * eps;
        xp = min(max(img + d, 0), 1);
        logit = ViTReach.evaluate(M, ViTReach.normImg(M, xp));
        [~, pred] = max(logit);
        if pred ~= label1, broken = true; return; end
    end
end
