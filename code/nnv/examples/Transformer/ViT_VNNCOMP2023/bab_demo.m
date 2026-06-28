function bab_demo()
% Show ReLU-phase-split BaB certifying where LP single-shot reach cannot.
here = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(here,'..','..','..','engine'))); addpath(here);
M = ViTReach.load(fullfile(here,'ibp_3_3_8.mat'));
k = 11; img = squeeze(M.images(k,:,:,:)); label = M.labels(k);
for eps255 = [0.80 0.90 1.00]
    eps = eps255/255;
    % single-shot LP (root, no splits)
    t0 = tic;
    ssLP = ViTReach.verifyBoxSplits(M, ViTReach_lb(M,img,eps), ViTReach_ub(M,img,eps), label, ...
                                    struct('mode','estimate','relu','fast'), struct('block',{},'neuron',{},'phase',{}));
    tss = toc(t0);
    % ReLU-split BaB
    t0 = tic;
    [bab, info] = ViTReach.verifyBaBRelu(M, img, label, eps, struct('babMaxNodes',31,'maxCand',12));
    tb = toc(t0);
    fprintf('eps=%.2f/255 | single-shot-LP=%d (%.1fs) | BaB=%d nodes=%d (%.1fs)\n', ...
        eps255, ssLP, tss, bab, info.nodes, tb);
end
end
function lb = ViTReach_lb(M,img,eps), [lb,~]=ViTReach.epsBox(M,img,eps); end
function ub = ViTReach_ub(M,img,eps), [~,ub]=ViTReach.epsBox(M,img,eps); end
