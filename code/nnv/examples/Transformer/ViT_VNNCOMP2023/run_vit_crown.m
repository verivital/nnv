function summary = run_vit_crown(modelName, nInst, opt)
%RUN_VIT_CROWN  Sound, LP-free verification of the VNN-COMP 2023 ViT instances.
%   Verifies argmax robustness at L-inf eps = 1/255 for every instance of one ViT
%   model using ViTCrown's native (from-scratch MATLAB; no external verifier) pipeline:
%       (1) lower the ViT to an op-DAG (ViTCrown.toOps),
%       (2) CROWN intermediate-bound refinement (ViTCrown.refineBounds), and
%       (3) alpha-optimization of the relu + attention-McCormick relaxations
%           (ViTCrown.optimizeAlpha).
%   No LP solver is used; every step is dense matrix algebra + gradient ascent.
%   The verdict is SOUND: a "robust" result means all 9 class-margin lower bounds
%   are > 0; otherwise the instance is reported "unknown" (the worst case).
%
%   Usage:
%       summary = run_vit_crown('ibp_3_3_8')             % all instances
%       summary = run_vit_crown('ibp_3_3_8', 20)         % first 20 instances
%       summary = run_vit_crown('pgd_2_3_16', inf, struct('nIter',150))
%
%   Options (opt struct):
%       eps        L-inf radius (default 1/255)
%       maxRefine  CROWN refinement passes (default 2)
%       nIter      alpha-opt iterations (default 100)
%       lr         alpha-opt learning rate (default 0.03)
%       selfcheck  random falsification self-check on every "robust" (default true)
%
%   The model bundle <modelName>.mat is produced by extract_weights.py from the
%   benchmark ONNX (gitignored, NOT committed). If missing, this errors with the
%   regeneration command. Depends only on ViTCrown.m + the (upstream) ViTReach and
%   SoftmaxAttn classes.

    if nargin < 2 || isempty(nInst), nInst = inf; end
    if nargin < 3, opt = struct(); end
    eps       = getf(opt,'eps',1/255);
    maxRefine = getf(opt,'maxRefine',2);
    nIter     = getf(opt,'nIter',100);
    lr        = getf(opt,'lr',0.03);
    selfcheck = getf(opt,'selfcheck',true);

    here = fileparts(mfilename('fullpath'));
    addpath(genpath(fullfile(here,'..','..','..','engine')));
    addpath(here);

    matPath = fullfile(here, [modelName '.mat']);
    assert(isfile(matPath), ['run_vit_crown:modelMissing  %s.mat not found.\n' ...
        'Regenerate the gitignored model bundle:  python extract_weights.py'], modelName);
    M   = ViTReach.load(matPath);
    ops = ViTCrown.toOps(M);                       % lower once (instance-independent)

    nImg  = size(M.images,1);
    nInst = min(nInst, nImg);
    margins = nan(nInst,1); times = nan(nInst,1); robust = false(nInst,1);

    fprintf('\n=== run_vit_crown %s : %d instances (LP-free, eps=%.5g) ===\n', modelName, nInst, eps);
    tAll = tic;
    for k = 1:nInst
        img = squeeze(M.images(k,:,:,:));          % [3,32,32] in [0,1]
        lab = M.labels(k);                         % 0-indexed
        t   = lab + 1; others = setdiff(1:10, t);
        C = zeros(9,10); for r = 1:9, C(r,t) = 1; C(r,others(r)) = -1; end
        [lb,ub] = ViTReach.epsBox(M, img, eps);

        tk = tic;
        [cl,cu,mm] = ViTCrown.refineBounds(ops, lb, ub, maxRefine);
        mg = ViTCrown.optimizeAlpha(ops, lb, ub, cl, cu, C, struct('nIter',nIter,'lr',lr));
        times(k) = toc(tk);
        margins(k) = min(mg);
        robust(k)  = all(mg > 0);

        % soundness self-check: a "robust" verdict must not be falsifiable
        if robust(k) && selfcheck
            if falsifiable(M, ops, img, lab, eps, C, 500)
                error('run_vit_crown:SOUNDNESS', ...
                    'instance %d certified robust but a counterexample was found!', k);
            end
        end
        fprintf('  inst %3d: label=%d  minMargin=%+.5f  robust=%d  (%.1fs)\n', ...
            k, lab, margins(k), robust(k), times(k));
    end
    nver = sum(robust);
    fprintf('--- run_vit_crown %s: %d/%d VERIFIED @ eps=%.5g; total %.0fs ---\n', ...
        modelName, nver, nInst, eps, toc(tAll));
    summary = struct('model',modelName,'nInst',nInst,'verified',nver, ...
        'robustMask',robust,'margins',margins,'times',times);
end

function tf = falsifiable(M, ops, img, lab, eps, C, nTry)
    tf = false;
    [lb,ub] = ViTReach.epsBox(M, img, eps);
    for i = 1:nTry
        xs = lb + (ub-lb).*rand(size(lb));         % sample the exact normalised eps-box
        if min(C * ViTCrown.evalOps(ops, xs)) < 0, tf = true; return; end
    end
end

function v = getf(s,f,d), if isfield(s,f), v = s.(f); else, v = d; end, end
