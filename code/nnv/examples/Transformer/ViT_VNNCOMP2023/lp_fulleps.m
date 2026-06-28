function lp_fulleps(modelName, nInst)
% Full-eps (1/255) verification using LP margins (marginMode='lp'), which honour
% the full Star constraint system (av-envelope facets + ReLU triangles) and are
% strictly tighter than the estimate margins used in run_benchmark. Reports the
% LP-verified count at full eps and the LP min-margin per instance.
    if nargin < 2, nInst = 15; end
    here = fileparts(mfilename('fullpath'));
    addpath(genpath(fullfile(here,'..','..','..','engine'))); addpath(here);
    M = ViTReach.load(fullfile(here, [modelName '.mat']));
    nInst = min(nInst, size(M.images,1));
    opt = struct('mode','estimate','relu','fast','marginMode','lp');
    verified = 0;
    fprintf('=== %s full-eps LP margins ===\n', modelName);
    for k = 1:nInst
        img = squeeze(M.images(k,:,:,:));
        [lb,ub] = ViTReach.epsBox(M, img, 1/255);
        t0 = tic;
        [rob, margins] = ViTReach.verify(M, lb, ub, M.labels(k), opt);
        verified = verified + (rob==1);
        fprintf('  inst %2d: LP robust=%d minMargin=%+.4f (%.1fs)\n', k, rob==1, min(margins), toc(t0));
    end
    fprintf('--- %s: %d/%d LP-verified @ full eps=1/255 ---\n', modelName, verified, nInst);
end
