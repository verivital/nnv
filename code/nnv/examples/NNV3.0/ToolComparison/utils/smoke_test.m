function smoke_test()
%SMOKE_TEST End-to-end sanity check of the ToolComparison wiring.
%
%   Six checks (~30 seconds total):
%     1. NNV is on the path and matlab2nnv works
%     2. AIVL functions resolve (verifyNetworkRobustness, estimateNetworkOutputBounds)
%     3. matlab2nnv handles a hand-built FC dlnetwork
%     4. NNV reach with relax-star-range produces a valid result
%     5. estimateNetworkOutputBounds runs and returns finite bounds
%     6. tool_utils helpers serialize and read back a result row
%
%   Optional 7th check (run only if a real ACAS Xu ONNX is available):
%     7. matlab2nnv on an ACAS Xu ONNX (proves the ScalingLayer engine fix)
%
%   On success, prints "SMOKE TEST PASSED" and exits 0.
%   On failure, prints which step failed and exits 1.

    fprintf("\n========== ToolComparison SMOKE TEST ==========\n");
    fail = @(step, msg) (assert(false, sprintf("[FAIL %s] %s", step, msg)));

    % --- 1. NNV path -----------------------------------------------------
    try
        v = which('NNVVERSION');
        assert(~isempty(v), "NNVVERSION not on path");
        fprintf("[1] NNV:                 %s\n", v);
    catch ME
        fail("1", sprintf("NNV path / install issue: %s", ME.message));
    end

    % --- 2. AIVL entry points -------------------------------------------
    %   AIVL is optional. If the Support Package isn't staged, we skip the
    %   MW-side checks (steps 2 and 5) and the rest of the smoke test still
    %   verifies the NNV-only path that 'tools',{'nnv'} runs use.
    %   Note: which() finds stub .m files shipped with MATLAB even when the
    %   support package isn't installed; we have to actually probe.
    haveAIVL = false;
    try
        netProbe = dlnetwork([featureInputLayer(2, "Name", "in"); ...
                              fullyConnectedLayer(2, "Name", "fc", ...
                                  "Weights", eye(2,'single'), "Bias", zeros(2,1,'single'))]);
        XLp = dlarray([0;0], "CB"); XUp = dlarray([1;1], "CB");
        [~, ~] = estimateNetworkOutputBounds(netProbe, XLp, XUp); %#ok<ASGLU>
        haveAIVL = true;
    catch
        haveAIVL = false;
    end
    if haveAIVL
        v1 = which('verifyNetworkRobustness');
        v2 = which('estimateNetworkOutputBounds');
        fprintf("[2] AIVL verifyNet:      %s\n", v1);
        fprintf("    AIVL estimateBounds: %s\n", v2);
    else
        fprintf("[2] AIVL not installed -- skipping MW-side checks (steps 2, 5).\n");
        fprintf("    To enable: stage atva26-aivl.tar.gz and run scripts/toolbox_install.m\n");
    end

    % --- 3. Build a tiny 3-layer FC network directly -------------------
    rng(0);
    W1 = randn(4, 5, 'single'); b1 = zeros(4, 1, 'single');
    W2 = randn(3, 4, 'single'); b2 = zeros(3, 1, 'single');
    layers = [
        featureInputLayer(5, "Name", "in")
        fullyConnectedLayer(4, "Name", "fc1", "Weights", W1, "Bias", b1)
        reluLayer("Name", "relu1")
        fullyConnectedLayer(3, "Name", "fc2", "Weights", W2, "Bias", b2)];
    netMW = dlnetwork(layers);
    fprintf("[3] hand-built dlnetwork: %d layers, 5->4->3\n", numel(netMW.Layers));

    try
        netNNV = matlab2nnv(netMW);
        fprintf("    matlab2nnv OK\n");
    catch ME
        fail("3", sprintf("matlab2nnv failed: %s", ME.message));
    end

    % --- 4. NNV reach (relax-star-range) on a small L_inf box -----------
    XLower = [-0.5; -0.5; -0.5; -0.5; -0.5];
    XUpper = [ 0.5;  0.5;  0.5;  0.5;  0.5];
    try
        IS = Star(XLower, XUpper);
        reachOpt = struct('reachOption', "single", 'numCores', 1, ...
                          'reachMethod', 'relax-star-range', 'relaxFactor', 0.5);
        t = tic;
        R = netNNV.reach(IS, reachOpt);
        tNNV = toc(t);
        assert(~isempty(R), "NNV reach returned empty result");
        fprintf("[4] NNV relax-star-range: %d output sets, %.2fs\n", numel(R), tNNV);
    catch ME
        fail("4", sprintf("NNV reach failed: %s", ME.message));
    end

    % --- 5. estimateNetworkOutputBounds ---------------------------------
    if haveAIVL
        try
            XL = dlarray(XLower, "CB");
            XU = dlarray(XUpper, "CB");
            t = tic;
            [yL, yU] = estimateNetworkOutputBounds(netMW, XL, XU);
            tMW = toc(t);
            yL = extractdata(yL); yU = extractdata(yU);
            assert(all(isfinite(yL)) && all(isfinite(yU)), "estimateNetworkOutputBounds returned non-finite values");
            assert(all(yU >= yL),                          "estimateNetworkOutputBounds bounds inverted");
            fprintf("[5] estimateBounds:      yL=%s yU=%s, %.2fs\n", ...
                    mat2str(yL(:)', 4), mat2str(yU(:)', 4), tMW);
        catch ME
            fprintf(2, "[5] estimateBounds FAILED:\n%s\n", getReport(ME, 'extended'));
            fail("5", sprintf("estimateNetworkOutputBounds failed: %s", ME.message));
        end
    else
        fprintf("[5] skipped (AIVL not installed)\n");
    end

    % --- 6. tool_utils round-trip ---------------------------------------
    try
        u = tool_utils();
        tmpFile = fullfile(tempdir, sprintf("smoke_%d.mat", randi(1e9)));
        cleanup = onCleanup(@() exist(tmpFile, 'file') && delete(tmpFile));
        row = u.new_row("nnv", "smoke", "instance0", "verified", 1.234, "relax-star-range-50", 30);
        u.append_to_mat(tmpFile, row);
        assert(u.has_instance(tmpFile, "nnv", "smoke", "instance0", "relax-star-range-50"), ...
               "has_instance returned false after append");
        loaded = u.load_results(tmpFile);
        assert(height(loaded) == 1 && loaded.status == "verified", "round-trip mismatch");
        fprintf("[6] tool_utils:          round-trip OK\n");
        delete(tmpFile);
    catch ME
        fail("6", sprintf("tool_utils round-trip failed: %s", ME.message));
    end

    % --- 7. (Optional) ACAS Xu ScalingLayer fix --------------------------
    %   If the CAV'23 ACAS Xu ONNX directory is reachable, verify that the
    %   ScalingLayer engine fix lets matlab2nnv import a real ACAS network.
    %   Skipped silently if the assets aren't installed.
    acas_dirs = { ...
        fullfile(fileparts(mfilename('fullpath')), '..', 'NNV2.0', 'Submission', ...
                 'CAV2023', 'NNV_vs_MATLAB', 'acas', 'onnx'), ...
        '/home/matlab/nnv/code/nnv/examples/NNV2.0/Submission/CAV2023/NNV_vs_MATLAB/acas/onnx'};
    acasOnnx = '';
    for k = 1:numel(acas_dirs)
        d = acas_dirs{k};
        if isfolder(d)
            f = dir(fullfile(d, '*.onnx'));
            if ~isempty(f), acasOnnx = fullfile(f(1).folder, f(1).name); break; end
        end
    end
    if isempty(acasOnnx)
        fprintf("[7] skipped (ACAS Xu ONNX assets not found)\n");
    else
        try
            netRaw = importONNXNetwork(acasOnnx, InputDataFormats="BCSS");
            nScale = sum(arrayfun(@(L) isa(L, 'nnet.cnn.layer.ScalingLayer'), netRaw.Layers));
            netNNV = matlab2nnv(netRaw); %#ok<NASGU>
            fprintf("[7] ACAS Xu ONNX import: ScalingLayer x%d, matlab2nnv OK\n", nScale);
        catch ME
            fail("7", sprintf("ACAS ONNX import failed (ScalingLayer fix): %s", ME.message));
        end
    end

    fprintf("\n========== SMOKE TEST PASSED ==========\n\n");
end
