%% test_ViTReach_reach - end-to-end sound reachability soundness + timing
% Propagates a per-pixel eps-box through the full ViT on Star sets and checks
% that EVERY sampled concrete logit vector lies inside the reach enclosure
% (estimate-box containment -- the estimate box outer-bounds the reachable set,
% so any sample escaping it is a definite soundness violation). Reports timing
% and verification margins. Each cell is self-contained (re-loads its model).

here = fileparts(mfilename('fullpath'));
nnvroot = fullfile(here, '..', '..', '..');
addpath(genpath(fullfile(nnvroot, 'engine')));
addpath(here);
rng(20260626);

%% pgd_2_3_16 reach is SOUND at full eps=1/255 + margins
here = pdir();
M = ViTReach.load(fullfile(here, 'pgd_2_3_16.mat'));
img = squeeze(M.images(1,:,:,:));
[lb, ub] = ViTReach.epsBox(M, img, 1/255);
t0 = tic;
L = ViTReach.reach(M, lb, ub, struct('mode','estimate','relu','fast'));
fprintf('pgd reach (full eps) time = %.1fs, nVar=%d\n', toc(t0), L.nVar);
[Llb, Lub] = L.estimateRanges();
viol = 0;
for t = 1:40
  d = (2*rand(size(img))-1) * (1/255);
  xp = min(max(img + d, 0), 1);
  logit = ViTReach.evaluate(M, ViTReach.normImg(M, xp));
  if any(logit < Llb - 1e-5) || any(logit > Lub + 1e-5), viol = viol + 1; end
end
fprintf('pgd soundness violations: %d/40\n', viol);
assert(viol == 0, 'pgd reach UNSOUND: %d violations', viol);
center_logit = ViTReach.evaluate(M, ViTReach.normImg(M, img));
assert(all(center_logit >= Llb - 1e-5) && all(center_logit <= Lub + 1e-5), 'center not in box');
[robust, margins] = ViTReach.verify(M, lb, ub, M.labels(1), struct('marginMode','estimate'));
fprintf('pgd full-eps verify: robust=%d, min margin=%.3f\n', robust, min(margins));

%% ibp_3_3_8 reach is SOUND at full eps (the verifiable model; depth 3, 17 tokens)
here = pdir();
M = ViTReach.load(fullfile(here, 'ibp_3_3_8.mat'));
img = squeeze(M.images(1,:,:,:));
[lb, ub] = ViTReach.epsBox(M, img, 1/255);
t0 = tic;
L = ViTReach.reach(M, lb, ub, struct('mode','estimate','relu','fast'));
fprintf('ibp reach (full eps) time = %.1fs, nVar=%d\n', toc(t0), L.nVar);
[Llb, Lub] = L.estimateRanges();
viol = 0;
for t = 1:25
  d = (2*rand(size(img))-1) * (1/255);
  xp = min(max(img + d, 0), 1);
  logit = ViTReach.evaluate(M, ViTReach.normImg(M, xp));
  if any(logit < Llb - 1e-5) || any(logit > Lub + 1e-5), viol = viol + 1; end
end
fprintf('ibp soundness violations: %d/25\n', viol);
assert(viol == 0, 'ibp reach UNSOUND: %d violations', viol);
[robust, margins] = ViTReach.verify(M, lb, ub, M.labels(1), struct('marginMode','estimate'));
fprintf('ibp full-eps verify: robust=%d, min margin=%.3f\n', robust, min(margins));

% ---- local helpers ----
function d = pdir()
  d = fileparts(mfilename('fullpath'));
end
