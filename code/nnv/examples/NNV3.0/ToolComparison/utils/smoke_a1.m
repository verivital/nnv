% smoke_a1.m -- End-to-end smoke for A1+A2+A5 inside nnv3.0:r2025b.
%
% (1) Verify ACAS Xu p3 NNV approx-star runs (1 net, 1 alg) — proves the
%     'aivl' tool dispatch and find_asset_dir() work post-rename.
% (2) Verify MNIST-ResNet AIVL deep-poly runs (1 point, 1 eps) —
%     proves verifyNetworkRobustness(net, XL, XU, ytrueIdx) works in R2025b
%     under the renamed schema (tool='aivl', algorithm='deep-poly').

addpath(genpath('/home/matlab/nnv/code/nnv'));
addpath(genpath('/home/matlab/nnv/code/nnv/examples/NNV3.0/ToolComparison'));

fprintf('\n========== Smoke A1: ACAS NNV ==========\n');
try
    run_acas_rl_tll('benchmarks',{'acas_p3'}, 'tools',{'nnv'}, ...
                    'algorithms',{'approx-star'}, 'numNets',1, 'timeout',30);
    fprintf('  ACAS NNV smoke OK\n');
catch ME
    fprintf(2, '  ACAS NNV smoke FAILED: %s\n', ME.message);
    fprintf(2, '%s\n', getReport(ME, 'extended'));
end

fprintf('\n========== Smoke A1: ACAS AIVL estimate-bounds ==========\n');
try
    run_acas_rl_tll('benchmarks',{'acas_p3'}, 'tools',{'aivl'}, ...
                    'algorithms',{'estimate-bounds'}, 'numNets',1, 'timeout',30);
    fprintf('  ACAS AIVL smoke OK\n');
catch ME
    fprintf(2, '  ACAS AIVL smoke FAILED: %s\n', ME.message);
end

fprintf('\n========== Smoke A1: MNIST-ResNet AIVL deep-poly ==========\n');
try
    run_mnist_resnet('models',{'mnist_resnet8'}, 'tools',{'aivl'}, ...
                     'algorithms',{'deep-poly'}, 'numPoints',1, ...
                     'epsilons',[1/255], 'timeout',60);
    fprintf('  MNIST-ResNet AIVL smoke OK\n');
catch ME
    fprintf(2, '  MNIST-ResNet AIVL smoke FAILED: %s\n', ME.message);
end

fprintf('\n========== Smoke complete ==========\n');
