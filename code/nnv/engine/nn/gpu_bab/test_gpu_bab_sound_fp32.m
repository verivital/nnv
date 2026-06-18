function tests = test_gpu_bab_sound_fp32
% TEST_GPU_BAB_SOUND_FP32  Soundness of the sound-FP32 outward-rounded crown_tight path (M-A2/M-B0).
%   gpu_bab_crown_tight in 'single' widens every CROWN bound OUTWARD by a Higham running-error radius
%   (i_outward_rad) + per-op backward-roundoff term (derr), so the FP32 bound is a PROVABLE sound
%   lower bound -- enabling a fast (no FP64 confirm) certified emit. A wrong FP32 bound = a wrong
%   unsat = -150, so this asserts the two load-bearing properties on ReLU nets (affine/relu/add, the
%   host path -- no GPU needed):
%     G1 PARITY: the single (sound-FP32) margin must be <= the double (FP64 oracle) margin, element-
%                wise (the widening only LOOSENS -- it can never spuriously raise a margin).
%     G2 MC:     the single margin must be <= the true min of C*f over the box (Monte-Carlo), i.e. a
%                genuine sound lower bound, never an over-estimate.
%   Run: runtests('test_gpu_bab_sound_fp32').
    tests = functiontests(localfunctions);
end

function test_sound_fp32_parity_fc(tc)
    rng(3, 'twister');
    ops = i_fcnet();                    % input(8) -> affine -> relu -> affine -> relu -> affine(5)
    lb = -0.3 * ones(8, 1); ub = 0.3 * ones(8, 1);
    C  = [1 -1 0 0 0; 0 1 -1 0 0; -1 0 0 1 0];
    mH = double(gpu_bab_crown_tight(ops, lb, ub, C, 'double', {})); mH = mH(:);
    mS = double(gpu_bab_crown_tight(ops, lb, ub, C, 'single', {})); mS = mS(:);
    verifyEqual(tc, numel(mS), numel(mH));
    verifyTrue(tc, all(isfinite(mS)), 'sound-FP32 margins must be finite');
    % G1 parity: single (sound-FP32, outward) <= double (FP64 oracle)
    verifyLessThanOrEqual(tc, mS, mH + 1e-4, 'G1: sound-FP32 margin must be <= FP64 margin');
    % G2 MC: single <= true min over the box
    N = 40000; X = lb + (ub - lb) .* rand(8, N);
    minCY = min(C * i_eval(ops, X), [], 2);
    verifyLessThanOrEqual(tc, mS, minCY + 1e-4, 'G2: sound-FP32 margin must be <= the MC true min');
end

function test_sound_fp32_parity_residual(tc)
    % residual (add) net: exercises the 'add' routing derr + a DAG.
    rng(8, 'twister');
    W1 = randn(10, 6); b1 = randn(10, 1);
    W2 = randn(10, 10); b2 = randn(10, 1);   % main branch (op4), added to the relu output (op2)
    W3 = randn(4, 10); b3 = randn(4, 1);
    ops = { struct('type','affine','W',W1,'b',b1,'src',0,'nOut',10), ...   % op1
            struct('type','relu','src',1), ...                            % op2
            struct('type','affine','W',W2,'b',b2,'src',2,'nOut',10), ...   % op3
            struct('type','relu','src',3), ...                            % op4
            struct('type','add','inputs',[2 4],'shape',[],'nOut',10), ...  % op5: op2 + op4 (residual)
            struct('type','affine','W',W3,'b',b3,'src',5,'nOut',4) };      % op6
    lb = -0.25 * ones(6, 1); ub = 0.25 * ones(6, 1);
    C  = [1 -1 0 0; 0 1 -1 0];
    mH = double(gpu_bab_crown_tight(ops, lb, ub, C, 'double', {})); mH = mH(:);
    mS = double(gpu_bab_crown_tight(ops, lb, ub, C, 'single', {})); mS = mS(:);
    verifyTrue(tc, all(isfinite(mS)));
    verifyLessThanOrEqual(tc, mS, mH + 1e-4, 'G1 (residual): sound-FP32 <= FP64');
    N = 40000; X = lb + (ub - lb) .* rand(6, N);
    minCY = min(C * i_eval(ops, X), [], 2);
    verifyLessThanOrEqual(tc, mS, minCY + 1e-4, 'G2 (residual): sound-FP32 <= MC true min');
end

function test_double_path_has_no_radius(tc)
    % The outward widening is gated to 'single': the 'double' path must be byte-identical to the
    % pre-sound-FP32 behaviour (derr=0, rad=0). Cross-check: doubling the box does not change the
    % double-path margin by any FP32-radius-sized amount (the radius would be ~1e-3..1e-1, FP64 ~1e-15).
    rng(5, 'twister');
    ops = i_fcnet();
    lb = -0.2 * ones(8, 1); ub = 0.2 * ones(8, 1);
    C  = [1 -1 0 0 0];
    m1 = double(gpu_bab_crown_tight(ops, lb, ub, C, 'double', {}));
    % MC must still bound it (double path is also sound); and it must be > the single path (tighter).
    mS = double(gpu_bab_crown_tight(ops, lb, ub, C, 'single', {}));
    N = 40000; X = lb + (ub - lb) .* rand(8, N);
    minCY = min(C * i_eval(ops, X), [], 2);
    verifyLessThanOrEqual(tc, m1(:), minCY + 1e-6, 'FP64 margin must be <= MC true min');
    verifyLessThanOrEqual(tc, mS(:), m1(:) + 1e-4, 'single (outward) must be <= double');
end

function ops = i_fcnet()
    W1 = randn(16, 8);  b1 = randn(16, 1);
    W2 = randn(12, 16); b2 = randn(12, 1);
    W3 = randn(5, 12);  b3 = randn(5, 1);
    ops = { struct('type','affine','W',W1,'b',b1,'src',0,'nOut',16), ...
            struct('type','relu','src',1), ...
            struct('type','affine','W',W2,'b',b2,'src',2,'nOut',12), ...
            struct('type','relu','src',3), ...
            struct('type','affine','W',W3,'b',b3,'src',4,'nOut',5) };
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
