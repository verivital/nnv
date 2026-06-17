function tests = test_gpu_bab_concat
% TEST_GPU_BAB_CONCAT  Soundness of the 'concat' op across the GPU-BaB engines.
%   Builds a small two-branch DAG (input -> two parallel affines -> concat -> relu -> affine),
%   then checks that gpu_bab_ibp / gpu_bab_crown_tight / gpu_bab_crown_spec_dag all produce
%   SOUND bounds (a guaranteed over-approximation) vs a dense Monte-Carlo ground truth, and
%   that the three engines agree at the root (no fixings). Concat is LINEAR (exact stacking),
%   so the only requirement is that the backward pass slices the coefficient back to the right
%   input block -- this test exercises exactly that (the second branch is a skip from the input,
%   forcing the full DAG path). Run: runtests('test_gpu_bab_concat').
    tests = functiontests(localfunctions);
end

function test_concat_sound(tc)
    rng(0, 'twister');
    n  = 4;
    lb = [-1; -0.5;  0; -2];
    ub = [ 1;  0.5;  1;  0.3];

    % branch 1: affine 3x4 ; branch 2: affine 2x4 (both from the INPUT, op 0) ;
    % concat -> 5 ; relu ; affine 2x5 -> out (2).
    A1 = randn(3, n); b1 = randn(3, 1);
    A2 = randn(2, n); b2 = randn(2, 1);
    A3 = randn(2, 5); b3 = randn(2, 1);
    ops = { struct('type','affine','W',A1,'b',b1,'src',0), ...
            struct('type','affine','W',A2,'b',b2,'src',0), ...
            struct('type','concat','inputs',[1 2],'sizes',[3 2],'src',1), ...
            struct('type','relu','src',3), ...
            struct('type','affine','W',A3,'b',b3,'src',4) };

    % spec: two linear output functionals (margins lower-bound C*out over the box)
    C = [1 -1; -1 1];

    % --- Monte-Carlo ground truth (exact forward eval) ---
    N = 40000;
    X = lb + (ub - lb) .* rand(n, N);
    Y = i_eval(ops, X);                       % nOut x N
    yMin = min(Y, [], 2); yMax = max(Y, [], 2);
    CYmin = min(C * Y, [], 2);

    tol = 1e-9;

    % --- IBP: sound interval over-approx of the output box ---
    [olb, oub] = gpu_bab_ibp(ops, lb, ub, 'double');
    verifyEqual(tc, numel(olb), numel(yMin), 'IBP output width mismatch (concat size wrong)');
    verifyLessThanOrEqual(tc, olb,            yMin + tol);   % lower bound is sound
    verifyGreaterThanOrEqual(tc, oub,         yMax - tol);   % upper bound is sound

    % --- crown_tight: sound lower bound on C*out ---
    m = gpu_bab_crown_tight(ops, lb, ub, C, 'double', {});
    m = m(:);
    verifyEqual(tc, numel(m), size(C,1), 'crown_tight margin count wrong');
    verifyTrue(tc, all(isfinite(m)), 'crown_tight produced non-finite margins');
    verifyLessThanOrEqual(tc, m,  CYmin + tol);              % margin <= true min (sound)

    % --- spec_dag (root, no fixings): sound AND consistent with crown_tight's IBP-forward path ---
    md = gpu_bab_crown_spec_dag(ops, lb, ub, C, 'double', {});
    md = md(:);
    verifyLessThanOrEqual(tc, md, CYmin + tol);              % sound

    % tightness sanity: bounds should not be absurdly loose (within, say, 50 of the MC range)
    verifyGreaterThan(tc, min(m), min(CYmin) - 50);
end

function Y = i_eval(ops, X)
% Exact forward evaluation of the op list (ground truth). Supports affine/relu/concat/add.
    cache = cell(numel(ops) + 1, 1);
    cache{1} = X;
    for k = 1:numel(ops)
        op = ops{k};
        switch op.type
            case 'affine'
                cache{k+1} = op.W * cache{op.src+1} + op.b(:);
            case 'relu'
                cache{k+1} = max(cache{op.src+1}, 0);
            case 'concat'
                cache{k+1} = cat(1, cache{op.inputs+1});
            case 'add'
                cache{k+1} = cache{op.inputs(1)+1} + cache{op.inputs(2)+1};
            otherwise
                error('i_eval:op', 'unsupported op "%s" in test ground-truth eval', op.type);
        end
    end
    Y = cache{end};
end
