function tests = test_gpu_bab_product
% TEST_GPU_BAB_PRODUCT  Soundness of the bilinear 'product' op (McCormick) across the GPU-BaB
%   engines. Builds input -> two parallel affines -> elementwise product -> affine, and checks
%   that gpu_bab_ibp / gpu_bab_crown_tight / gpu_bab_crown_spec_dag produce SOUND bounds vs a
%   dense Monte-Carlo ground truth of the (nonlinear) network. Also directly verifies that the
%   gpu_bab_mul_relax planes envelope x.*y over a box. The product is the genuine nonlinearity
%   in lsnc Lyapunov nets (x .* (P*x)); McCormick has an irreducible gap (closed only by
%   branching), so this checks SOUNDNESS, not tightness. Run: runtests('test_gpu_bab_product').
    tests = functiontests(localfunctions);
end

function test_product_sound(tc)
    rng(1, 'twister');
    n  = 4;
    lb = [-1; -0.5; -0.8; -1.2];
    ub = [ 1;  0.5;  0.8;  0.7];
    A1 = randn(3, n); b1 = randn(3, 1);
    A2 = randn(3, n); b2 = randn(3, 1);
    A3 = randn(2, 3); b3 = randn(2, 1);
    ops = { struct('type','affine','W',A1,'b',b1,'src',0), ...
            struct('type','affine','W',A2,'b',b2,'src',0), ...
            struct('type','product','inputs',[1 2],'sizes',[3 3],'src',1), ...
            struct('type','affine','W',A3,'b',b3,'src',3) };
    C = [1 -1; -1 1];
    tol = 1e-9;

    % --- Monte-Carlo ground truth (exact nonlinear forward eval) ---
    N = 60000;
    X = lb + (ub - lb) .* rand(n, N);
    Y = i_eval(ops, X);
    yMin = min(Y, [], 2); yMax = max(Y, [], 2);
    CYmin = min(C * Y, [], 2);

    % --- IBP: sound interval over the nonlinear output ---
    [olb, oub] = gpu_bab_ibp(ops, lb, ub, 'double');
    verifyEqual(tc, numel(olb), numel(yMin), 'IBP output width mismatch');
    verifyLessThanOrEqual(tc, olb, yMin + tol);
    verifyGreaterThanOrEqual(tc, oub, yMax - tol);

    % --- crown_tight: sound lower bound on C*out (McCormick relaxed) ---
    m = gpu_bab_crown_tight(ops, lb, ub, C, 'double', {}); m = m(:);
    verifyEqual(tc, numel(m), size(C,1), 'crown_tight margin count wrong');
    verifyTrue(tc, all(isfinite(m)), 'crown_tight produced non-finite margins');
    verifyLessThanOrEqual(tc, m, CYmin + tol);

    % --- spec_dag (root, no fixings): sound ---
    md = gpu_bab_crown_spec_dag(ops, lb, ub, C, 'double', {}); md = md(:);
    verifyLessThanOrEqual(tc, md, CYmin + tol);

    % --- direct McCormick-envelope soundness (planes bracket x.*y over the box) ---
    xl = [-2; -1; 0.5]; xu = [1; 3; 2];
    yl = [-1; -2; -0.5]; yu = [2; 1; 1.5];
    [aL,bL,cL,aU,bU,cU] = gpu_bab_mul_relax(xl, xu, yl, yu, [], [], 'double');
    M = 8000;
    for i = 1:numel(xl)
        xs = xl(i) + (xu(i) - xl(i)) * rand(1, M);
        ys = yl(i) + (yu(i) - yl(i)) * rand(1, M);
        z  = xs .* ys;
        lo = aL(i)*xs + bL(i)*ys + cL(i);
        hi = aU(i)*xs + bU(i)*ys + cU(i);
        verifyLessThanOrEqual(tc, lo, z + tol, sprintf('lower plane not an under-estimator (elt %d)', i));
        verifyGreaterThanOrEqual(tc, hi, z - tol, sprintf('upper plane not an over-estimator (elt %d)', i));
    end
end

function Y = i_eval(ops, X)
% Exact forward evaluation (ground truth). Supports affine/relu/product/concat/add.
    cache = cell(numel(ops) + 1, 1); cache{1} = X;
    for k = 1:numel(ops)
        op = ops{k};
        switch op.type
            case 'affine',  cache{k+1} = op.W * cache{op.src+1} + op.b(:);
            case 'relu',    cache{k+1} = max(cache{op.src+1}, 0);
            case 'product', cache{k+1} = cache{op.inputs(1)+1} .* cache{op.inputs(2)+1};
            case 'concat',  cache{k+1} = cat(1, cache{op.inputs+1});
            case 'add',     cache{k+1} = cache{op.inputs(1)+1} + cache{op.inputs(2)+1};
            otherwise, error('i_eval:op', 'unsupported op "%s"', op.type);
        end
    end
    Y = cache{end};
end
