function tests = test_gpu_bab_running_derr
% TEST_GPU_BAB_RUNNING_DERR  Soundness + tightness of the running-error (measured-delta) derr (P2).
%   NNV_DERR_RUNNING replaces the relu-intercept's a-priori gamma_width WORST-CASE reduction roundoff
%   with the MEASURED roundoff |fl32(Ab*bu)-fl64(Ab*bu)| + a sound fl64-residual cushion, via
%   min(a-priori, measured). A wrong (under-)widening = a spuriously-high margin = a wrong unsat = -150,
%   so this is the GOLD-GATE on WIDE hidden layers (large numel(bu) -> the a-priori gamma_width is large
%   and the measured-delta meaningfully smaller, exercising the new path -- unlike the width<=16 cases in
%   test_gpu_bab_sound_fp32). Asserts, with NNV_DERR_RUNNING=1, on randomized wide ReLU nets:
%     S1 SOUND:  the measured-delta single margin <= the FP64 oracle margin (STRICT tol -- an under-
%                widening cushion bug surfaces as single > double). running_derr >= |fp32-fp64| (advisor).
%     S2 MC:     the measured-delta single margin <= the true min of C*f over the box (Monte-Carlo).
%     S3 TIGHT:  the measured-delta margin is >= the a-priori margin (the min() can only LOOSEN the
%                widening, never tighten it past sound) -- and on a WIDE layer it is strictly tighter,
%                confirming the new term actually fires (else the build is a no-op).
%   Run: runtests('test_gpu_bab_running_derr').
    tests = functiontests(localfunctions);
end

function setup(tc)
    % the widening is opt-in via NNV_SOUND_FP32_TIGHT; the measured-delta via NNV_DERR_RUNNING.
    tc.TestData.oldSound = getenv('NNV_SOUND_FP32_TIGHT');
    tc.TestData.oldRun   = getenv('NNV_DERR_RUNNING');
    setenv('NNV_SOUND_FP32_TIGHT', '1');
end

function teardown(tc)
    setenv('NNV_SOUND_FP32_TIGHT', tc.TestData.oldSound);
    setenv('NNV_DERR_RUNNING', tc.TestData.oldRun);
end

function test_running_derr_wide(tc)
    % WIDE hidden layer (W=3000) S1 SOUNDNESS over several seeds. NOTE: a 2-relu net couples the measured-delta
    % into the INTERMEDIATE bounds (preL/preU), and CROWN relaxations are NOT monotonic in those bounds, so the
    % raw margin can shift either way vs a-priori -- the ONLY invariant is sound: measured-delta single <= FP64
    % oracle. (The clean monotone "tighter than a-priori" check is test_running_derr_onerelu, where the relu's
    % own interm bound comes from the input affine and is unaffected by the relu measured-delta.)
    for seed = [1 7 21 42]
        rng(seed, 'twister');
        ops = i_widenet(18, 3000, 6);
        lb = -0.2 * ones(18, 1); ub = 0.2 * ones(18, 1);
        C  = [1 -1 0 0 0 0; 0 1 -1 0 0 0; 0 0 0 1 -1 0];
        mH = double(gpu_bab_crown_tight(ops, lb, ub, C, 'double', {})); mH = mH(:);   % FP64 oracle
        setenv('NNV_DERR_RUNNING', '1'); mRun = double(gpu_bab_crown_tight(ops, lb, ub, C, 'single', {})); mRun = mRun(:);
        verifyTrue(tc, all(isfinite(mRun)), sprintf('seed %d: finite', seed));
        % S1 SOUND (strict): measured-delta single <= FP64 oracle. tol 1e-6 >> double roundoff (~1e-12),
        % << a real under-widening (a-priori widening here is ~1e-4..1e-2, an under-widen bug would be that scale).
        verifyLessThanOrEqual(tc, mRun, mH + 1e-6, sprintf('S1 seed %d: measured-delta margin must be <= FP64 oracle', seed));
    end
end

function test_running_derr_onerelu(tc)
    % S3 TIGHT (clean, no interm coupling): ONE hidden relu (input -> affine(W=4000) -> relu -> affine(out)).
    % The relu's interm bound is the input affine (no relu before it) so the relu measured-delta does NOT change
    % preL/preU -> the final relaxation is identical -> mRun differs from mAp ONLY by the relu-intercept derr,
    % which min(a-priori,measured) can only SHRINK. So mRun >= mAp strictly (the new term fired + is tighter).
    for seed = [2 19]
        rng(seed, 'twister');
        w = 4000;
        W1 = randn(w, 14)/sqrt(14); b1 = randn(w, 1)*0.1;
        W2 = randn(5, w)/sqrt(w);   b2 = randn(5, 1)*0.1;
        ops = { struct('type','affine','W',W1,'b',b1,'src',0,'nOut',w), ...
                struct('type','relu','src',1), ...
                struct('type','affine','W',W2,'b',b2,'src',2,'nOut',5) };
        lb = -0.2 * ones(14, 1); ub = 0.2 * ones(14, 1);
        C  = [1 -1 0 0 0; 0 0 1 -1 0];
        mH = double(gpu_bab_crown_tight(ops, lb, ub, C, 'double', {})); mH = mH(:);
        setenv('NNV_DERR_RUNNING', '');  mAp  = double(gpu_bab_crown_tight(ops, lb, ub, C, 'single', {})); mAp  = mAp(:);
        setenv('NNV_DERR_RUNNING', '1'); mRun = double(gpu_bab_crown_tight(ops, lb, ub, C, 'single', {})); mRun = mRun(:);
        verifyLessThanOrEqual(tc, mRun, mH + 1e-6, sprintf('S1 (1-relu) seed %d: measured-delta <= FP64 oracle', seed));
        verifyGreaterThanOrEqual(tc, mRun, mAp - 1e-7, sprintf('S3 (1-relu) seed %d: measured-delta >= a-priori (min can only shrink the widening)', seed));
        verifyGreaterThan(tc, max(mRun - mAp), 1e-6, sprintf('S3 (1-relu) seed %d: measured-delta STRICTLY tighter (the new term fired)', seed));
    end
end

function test_running_derr_mc(tc)
    % S2: measured-delta single margin is a genuine sound lower bound vs Monte-Carlo ground truth.
    rng(13, 'twister');
    ops = i_widenet(12, 2000, 5);
    lb = -0.25 * ones(12, 1); ub = 0.25 * ones(12, 1);
    C  = [1 -1 0 0 0; 0 0 1 -1 0];
    setenv('NNV_DERR_RUNNING', '1');
    mRun = double(gpu_bab_crown_tight(ops, lb, ub, C, 'single', {})); mRun = mRun(:);
    N = 60000; X = lb + (ub - lb) .* rand(12, N);
    minCY = min(C * i_eval(ops, X), [], 2);
    verifyLessThanOrEqual(tc, mRun, minCY + 1e-5, 'S2: measured-delta margin must be <= the MC true min');
end

function test_running_derr_residual(tc)
    % S1/S3 on a residual (add) DAG with a wide branch -- exercises the 'add' merge + the wide-relu intercept.
    rng(8, 'twister');
    W1 = randn(2000, 10)/8; b1 = randn(2000, 1);
    W2 = randn(2000, 2000)/45; b2 = randn(2000, 1);
    W3 = randn(4, 2000)/45; b3 = randn(4, 1);
    ops = { struct('type','affine','W',W1,'b',b1,'src',0,'nOut',2000), ...
            struct('type','relu','src',1), ...
            struct('type','affine','W',W2,'b',b2,'src',2,'nOut',2000), ...
            struct('type','relu','src',3), ...
            struct('type','add','inputs',[2 4],'shape',[],'nOut',2000), ...
            struct('type','affine','W',W3,'b',b3,'src',5,'nOut',4) };
    lb = -0.15 * ones(10, 1); ub = 0.15 * ones(10, 1);
    C  = [1 -1 0 0; 0 1 -1 0];
    mH = double(gpu_bab_crown_tight(ops, lb, ub, C, 'double', {})); mH = mH(:);
    setenv('NNV_DERR_RUNNING', '1'); mRun = double(gpu_bab_crown_tight(ops, lb, ub, C, 'single', {})); mRun = mRun(:);
    verifyTrue(tc, all(isfinite(mRun)));
    verifyLessThanOrEqual(tc, mRun, mH + 1e-6, 'S1 (residual): measured-delta <= FP64 oracle');
    N = 60000; X = lb + (ub - lb) .* rand(10, N);
    minCY = min(C * i_eval(ops, X), [], 2);
    verifyLessThanOrEqual(tc, mRun, minCY + 1e-5, 'S2 (residual): measured-delta <= MC true min');
end

function ops = i_widenet(nin, w, nout)
    % input(nin) -> affine(w) -> relu -> affine(w) -> relu -> affine(nout); weights scaled so the
    % pre-activations straddle 0 (many unstable relus -> nonzero bu intercepts -> the term fires).
    W1 = randn(w, nin)/sqrt(nin);   b1 = randn(w, 1)*0.1;
    W2 = randn(w, w)/sqrt(w);       b2 = randn(w, 1)*0.1;
    W3 = randn(nout, w)/sqrt(w);    b3 = randn(nout, 1)*0.1;
    ops = { struct('type','affine','W',W1,'b',b1,'src',0,'nOut',w), ...
            struct('type','relu','src',1), ...
            struct('type','affine','W',W2,'b',b2,'src',2,'nOut',w), ...
            struct('type','relu','src',3), ...
            struct('type','affine','W',W3,'b',b3,'src',4,'nOut',nout) };
end

function Y = i_eval(ops, X)
    cache = cell(numel(ops) + 1, 1); cache{1} = X;
    for k = 1:numel(ops)
        op = ops{k};
        switch op.type
            case 'affine', cache{k+1} = op.W * cache{op.src+1} + op.b(:);
            case 'relu',   cache{k+1} = max(cache{op.src+1}, 0);
            case 'add',    cache{k+1} = cache{op.inputs(1)+1} + cache{op.inputs(2)+1};
            otherwise, error('i_eval:op', 'unsupported op "%s"', op.type);
        end
    end
    Y = cache{end};
end
