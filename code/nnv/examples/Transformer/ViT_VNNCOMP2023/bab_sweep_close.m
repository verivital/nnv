function bab_sweep_close(K, budget)
% Efficient full-eps BaB sweep on the instances that could plausibly verify: rank
% all 15 ibp_3_3_8 instances by estimate margin (cheap), then run FF-ReLU-split BaB
% at eps=1/255 on the closest K. Reports, per instance, the single-shot LP margin
% (root, honours the av-envelope facets) and the BaB verdict, plus the total count.
% This is the practical "only run the verifiable subset" -- the closest-gap
% instances are exactly where alpha,beta-CROWN's 79 live; my weaker (attention
% box-lifting) bounder will not match 79, but this measures what it CAN do.
    if nargin < 1, K = 5; end
    if nargin < 2, budget = 40; end
    here = fileparts(mfilename('fullpath'));
    addpath(genpath(fullfile(here,'..','..','..','engine'))); addpath(here);
    M = ViTReach.load(fullfile(here,'ibp_3_3_8.mat'));
    n = 15; eps = 1/255;
    optE = struct('mode','estimate','relu','fast','marginMode','estimate');

    % rank by estimate margin (closest to 0 first)
    estM = nan(n,1);
    for k = 1:n
        img = squeeze(M.images(k,:,:,:));
        [lb,ub] = ViTReach.epsBox(M, img, eps);
        [~, margins] = ViTReach.verify(M, lb, ub, M.labels(k), optE);
        estM(k) = min(margins);
    end
    [~, ord] = sort(estM, 'descend');
    cand = ord(1:K)';
    fprintf('Closest %d ibp instances by estimate margin: %s\n', K, mat2str(cand));

    emptySplits = struct('block',{},'neuron',{},'phase',{});
    verified = 0;
    for k = cand
        img = squeeze(M.images(k,:,:,:)); label = M.labels(k);
        [lb,ub] = ViTReach.epsBox(M, img, eps);
        % single-shot LP (root) margin
        [~, rootLP] = ViTReach.verifyBoxSplits(M, lb, ub, label, optE, emptySplits);
        % full-eps BaB
        t0 = tic;
        [st, info] = ViTReach.verifyBaBRelu(M, img, label, eps, struct('babMaxNodes',budget,'maxCand',18));
        verified = verified + (st==1);
        fprintf('  inst %2d: estM=%+.4f rootLP=%+.4f | BaB status=%d nodes=%d (%.0fs)%s\n', ...
            k, estM(k), rootLP, st, info.nodes, toc(t0), tern(st==1,'  <-- VERIFIED @ full eps', ''));
    end
    fprintf('=== BaB @ full eps=1/255: %d/%d verified (closest %d ibp instances) ===\n', verified, K, K);
    fprintf('BAB_SWEEP_DONE\n');
end
function s = tern(c,a,b), if c, s=a; else, s=b; end, end
