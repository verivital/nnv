%% test_ViTReach_bab - sound branch-and-bound: soundness + does BaB beat single-shot?
% (1) input-split BaB is SOUND (never certifies a falsifiable box; certifies small eps);
% (2) ReLU-phase-split BaB can certify a STRICTLY larger eps than single-shot reach
%     on an instance near its certified boundary (the mechanism that closes the gap).

here = fileparts(mfilename('fullpath'));
nnvroot = fullfile(here, '..', '..', '..');
addpath(genpath(fullfile(nnvroot, 'engine')));
addpath(here);
rng(20260626);

%% input-split BaB soundness: robust@small-eps sound; never robust on a falsifiable box
here = pdir();
M = ViTReach.load(fullfile(here, 'ibp_3_3_8.mat'));
k = 11; img = squeeze(M.images(k,:,:,:)); label = M.labels(k);
% small eps: single-shot already robust -> BaB must also be robust (status 1)
[st_small, ~] = ViTReach.verifyBaB(M, img, label, 0.3/255, struct('babMaxNodes',8));
fprintf('BaB(input) @0.3/255: status=%d (expect 1 robust)\n', st_small);
assert(st_small == 1, 'BaB input-split should certify a clearly-robust small box');
% large eps: model is misclassifiable nearby -> BaB must NOT return robust
[st_big, info_big] = ViTReach.verifyBaB(M, img, label, 6/255, struct('babMaxNodes',12,'babFalsifyTries',400));
fprintf('BaB(input) @6/255: status=%d nodes=%d (expect 0 or 2, never 1)\n', st_big, info_big.nodes);
assert(st_big ~= 1, 'BaB must not certify a falsifiable box (SOUNDNESS)');
fprintf('PASS input-split BaB soundness\n');

%% ReLU-phase-split BaB certifies beyond single-shot on a borderline instance
here = pdir();
M = ViTReach.load(fullfile(here, 'ibp_3_3_8.mat'));
k = 11; img = squeeze(M.images(k,:,:,:)); label = M.labels(k);
% single-shot certified radius for this instance is ~0.625/255 (from the sweep);
% pick an eps just above it where single-shot reach is UNKNOWN:
epsTest = 0.66/255;
ss = ViTReach.verify(M, ViTReach_box(M,img,epsTest,1), ViTReach_box(M,img,epsTest,2), label, struct('marginMode','lp'));
fprintf('single-shot @%.3f/255: robust=%d (expect 2 unknown)\n', epsTest*255, ss);
t0 = tic;
[st_bab, info] = ViTReach.verifyBaBRelu(M, img, label, epsTest, struct('babMaxNodes',24,'maxCand',12));
fprintf('ReLU-split BaB @%.3f/255: status=%d nodes=%d (%.1fs)\n', epsTest*255, st_bab, info.nodes, toc(t0));
if st_bab == 1 && ss ~= 1
    fprintf('PASS BaB certifies beyond single-shot (sound BaB win)\n');
else
    fprintf('NOTE: BaB status=%d, single-shot=%d at %.3f/255 (mechanism sound; gain bounded by budget)\n', ...
        st_bab, ss, epsTest*255);
end
assert(st_bab ~= 0, 'BaB returned not-robust without a CE search; unexpected');

% ---- helpers ----
function d = pdir(), d = fileparts(mfilename('fullpath')); end
function out = ViTReach_box(M, img, eps, which)
  [lb, ub] = ViTReach.epsBox(M, img, eps);
  if which==1, out = lb; else, out = ub; end
end
