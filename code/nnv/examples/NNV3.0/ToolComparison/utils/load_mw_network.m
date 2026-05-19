function net = load_mw_network(onnx, inputFmt, outputFmt)
%LOAD_MW_NETWORK  R2025b-compatible ONNX loader.
%   Prefer importNetworkFromONNX (returns dlnetwork with affine folds).
%   Fall back to the manual ElementwiseAffineLayer fold for older opsets.
%   Accepts inputFmt in {'BCSS','BC'}. Optional outputFmt in {'BC',''}.
%
%   importNetworkFromONNX writes auto-generated custom layers to a
%   `+pkgname` directory in the current working directory. When the cwd
%   is read-only (e.g. a bind-mounted host source tree), the import
%   fails. We cd to tempdir for the import and restore the original cwd.
%
%   Copied from ToolComparison/acas_rl_tll/run_acas_rl_tll.m so v2 has
%   it as a top-level callable. Single source of truth lives in v2/utils/.

    if nargin < 3, outputFmt = ''; end
    oldDir   = pwd;
    importDir = tempname; mkdir(importDir);
    cleanupCd  = onCleanup(@() safe_cd(oldDir));         %#ok<NASGU>
    cleanupRm  = onCleanup(@() safe_rmdir(importDir));   %#ok<NASGU>
    cd(importDir);
    try
        if isempty(outputFmt)
            net = importNetworkFromONNX(onnx, InputDataFormats=string(inputFmt));
        else
            net = importNetworkFromONNX(onnx, ...
                InputDataFormats=string(inputFmt), OutputDataFormats=string(outputFmt));
        end
        return;
    catch
    end
    try
        if isempty(outputFmt)
            net = importONNXNetwork(onnx, InputDataFormats=string(inputFmt));
        else
            net = importONNXNetwork(onnx, ...
                InputDataFormats=string(inputFmt), OutputDataFormats=string(outputFmt));
        end
        Layers = net.Layers;
        n = numel(Layers);
        keep = [];
        for i = 1:(n-1)
            if isa(Layers(i), 'nnet.onnx.layer.ElementwiseAffineLayer') && i > 1
                Layers(i-1).Bias = Layers(i).Offset;
            else
                keep(end+1) = i; %#ok<AGROW>
            end
        end
        Layers = Layers(keep);
        net    = dlnetwork(Layers);
    catch ME2
        error('load_mw_network:import_failed', ...
              'Both importNetworkFromONNX and importONNXNetwork failed: %s', ME2.message);
    end
end

function safe_cd(d)
    if isfolder(d), cd(d); end
end

function safe_rmdir(d)
    if isfolder(d)
        try, rmdir(d, 's'); catch, end %#ok<NOCOM>
    end
end
