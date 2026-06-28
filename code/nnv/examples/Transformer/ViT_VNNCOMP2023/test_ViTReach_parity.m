%% test_ViTReach_parity - the MATLAB ViT forward must match onnxruntime
% Validates evaluate() (patch-embed layout, BN folding, multi-head attention,
% FF ReLU, meanpool, head) against the onnxruntime reference logits stored in the
% .mat by extract_weights.py. If this passes, reach() verifies the right function.

here = fileparts(mfilename('fullpath'));
nnvroot = fullfile(here, '..', '..', '..');
addpath(genpath(fullfile(nnvroot, 'engine')));
addpath(here);

%% pgd_2_3_16 parity vs onnxruntime
here = pdir();
M = ViTReach.load(fullfile(here, 'pgd_2_3_16.mat'));
ref = M.ref_logits;                       % [5 x 10]
maxerr = 0;
for k = 1:5
  img = squeeze(M.ref_norm_images(k,:,:,:));     % [3,32,32] normalised
  logits = ViTReach.evaluate(M, img);
  maxerr = max(maxerr, max(abs(logits(:) - ref(k,:)')));
end
fprintf('pgd_2_3_16 max parity err = %.3e\n', maxerr);
assert(maxerr < 1e-4, 'pgd parity FAILED: %.3e', maxerr);

%% ibp_3_3_8 parity vs onnxruntime
here = pdir();
M = ViTReach.load(fullfile(here, 'ibp_3_3_8.mat'));
ref = M.ref_logits;
maxerr = 0;
for k = 1:5
  img = squeeze(M.ref_norm_images(k,:,:,:));
  logits = ViTReach.evaluate(M, img);
  maxerr = max(maxerr, max(abs(logits(:) - ref(k,:)')));
end
fprintf('ibp_3_3_8 max parity err = %.3e\n', maxerr);
assert(maxerr < 1e-4, 'ibp parity FAILED: %.3e', maxerr);

% ---- local helpers ----
function d = pdir()
  d = fileparts(mfilename('fullpath'));
end
