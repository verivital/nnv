% probe_aivl_vnnlib.m -- Investigate what AIVL R2025b actually accepts.
%
% Phase A1 gate: we already know `verifyNetworkRobustness(net, vnnlibFile)`
% does NOT exist in R2025b (the VNNLIB ingest landed in R2026a). Probe the
% argmax-form signature on an ACAS Xu net: ACAS p3/p4 are expressible as
% argmax (output 1 should dominate), so a proper AIVL DeepPoly invocation
% is `verifyNetworkRobustness(net, XL, XU, classes=1)`.

addpath(genpath('/home/matlab/nnv/code/nnv'));
addpath(genpath('/home/matlab/nnv/code/nnv/examples/NNV3.0/ToolComparison'));

fprintf('=== AIVL on path (clean — aivnv dir only, no contents) ===\n');
matches = dir(fullfile(userpath(), 'SupportPackages', 'R*', 'toolbox', 'nnet', 'supportpackages', 'aivnv'));
if isempty(matches)
    fprintf(2, '  AIVL not installed; aborting probe.\n'); return;
end
addpath(matches(1).folder, '-begin');
fprintf('  added: %s\n', matches(1).folder);
fprintf('  verifyNetworkRobustness    -> %s\n', which('verifyNetworkRobustness'));
fprintf('  estimateNetworkOutputBounds -> %s\n', which('estimateNetworkOutputBounds'));

fprintf('\n=== Build ACAS Xu network ===\n');
acasDir = '/home/matlab/nnv/code/nnv/examples/NNV3.0/ToolComparison/acas_rl_tll/acas';
onnxFile = fullfile(acasDir, 'onnx', 'ACASXU_run2a_1_1_batch_2000.onnx');

oldDir = pwd;
tmpDir = tempname; mkdir(tmpDir); cd(tmpDir);
cleanup1 = onCleanup(@() cd(oldDir));
cleanup2 = onCleanup(@() rmdir(tmpDir, 's'));

netRaw = importNetworkFromONNX(onnxFile, InputDataFormats="BCSS");
netMW  = rebuild_for_aivl(netRaw);
fprintf('  network rebuilt: %d layers (FC stack)\n', numel(netMW.Layers));

fprintf('\n=== Test argmax-form verifyNetworkRobustness on ACAS p3 ===\n');
% ACAS Xu property 3: for the input box below, output 1 (COC) should NOT be
% the minimum -- equivalently, output 1 should dominate the others. Argmax
% form: `class=1` means "verify that argmax(net(x)) == 1 for all x in box".
% (NOTE: ACAS p3 actually says "output 1 is the minimum is a violation".
% We probe the API shape; semantic correctness is a Phase A1 follow-up.)
XLower = [-0.303531156; -0.009549297; 0.493380324; 0.3; 0.3];
XUpper = [-0.298552812;  0.009549297; 0.5;        0.5; 0.5];
XL = dlarray(single(XLower), "CB");
XU = dlarray(single(XUpper), "CB");

fprintf('\n  -- (1a) verifyNetworkRobustness(net, XL, XU, 1) [class=1, default alg]\n');
try
    t = tic;
    r = verifyNetworkRobustness(netMW, XL, XU, 1);
    fprintf('     OK in %.2fs, result class=%s, value=%s\n', toc(t), class(r), string(r));
catch ME
    fprintf(2, '     FAIL: %s\n', ME.message);
    fprintf(2, '     (id=%s)\n', ME.identifier);
end

fprintf('\n  -- (1b) verifyNetworkRobustness(net, XL, XU, 1, Algorithm="deep-poly")\n');
try
    t = tic;
    r = verifyNetworkRobustness(netMW, XL, XU, 1, Algorithm="deep-poly");
    fprintf('     OK in %.2fs, result=%s\n', toc(t), string(r));
catch ME
    fprintf(2, '     FAIL: %s\n', ME.message);
end

fprintf('\n  -- (1c) probe Algorithm name space\n');
try
    verifyNetworkRobustness(netMW, XL, XU, 1, Algorithm="__probe__");
catch ME
    fprintf('     probe error (lists allowed values): %s\n', ME.message);
end

fprintf('\n  -- (1d) Help signature\n');
try
    h = help('verifyNetworkRobustness');
    fprintf('%s\n', strtrim(strjoin(strsplit(h, newline), [newline '     '])));
catch
end

fprintf('\n=== Probe done ===\n');
