% smoke_a2.m -- end-to-end smoke for the A2 schema rename + A5 self-contained
% assets, running inside the nnv3.0:r2025b container.

addpath(genpath('/home/matlab/nnv/code/nnv'));
addpath(genpath('/home/matlab/nnv/code/nnv/examples/NNV3.0/ToolComparison'));

fprintf('=== canonicalize ===\n');
canonicalize_bundled_results();

fprintf('\n=== make_acas_rl_tll_table ===\n');
try
    make_acas_rl_tll_table();
    fprintf('  table_A OK\n');
catch ME
    fprintf(2, '  FAILED: %s\n', ME.message);
end

fprintf('\n=== make_mnist_resnet_table ===\n');
try
    make_mnist_resnet_table();
    fprintf('  table_C OK\n');
catch ME
    fprintf(2, '  FAILED: %s\n', ME.message);
end

fprintf('\n=== Asset locations check ===\n');
asset_subs = {'acas','rl_benchmarks','oval21'};
for k = 1:numel(asset_subs)
    sub = asset_subs{k};
    local = fullfile('/home/matlab/nnv/code/nnv/examples/NNV3.0/ToolComparison/acas_rl_tll', sub);
    if isfolder(local)
        n_onnx = numel(dir(fullfile(local,'onnx','*.onnx')));
        n_vnn  = numel(dir(fullfile(local,'vnnlib','*.vnnlib')));
        fprintf('  %s: %d onnx, %d vnnlib\n', sub, n_onnx, n_vnn);
    else
        fprintf(2, '  %s: NOT FOUND at %s\n', sub, local);
    end
end

fprintf('\n=== Smoke OK ===\n');
