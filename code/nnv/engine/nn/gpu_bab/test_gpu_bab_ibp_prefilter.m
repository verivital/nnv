function tests = test_gpu_bab_ibp_prefilter
% TEST_GPU_BAB_IBP_PREFILTER  Soundness of the IBP-prefilter for the sound-path interm bounds (perf lever).
%   NNV_CROWN_IBP_PREFILTER skips the chunked CROWN backward for relu neurons the SOUND DOUBLE IBP proves
%   stable (ibpLo>=0 or ibpHi<=0) and fills them from the (sound) IBP bound. A stable relu's relaxation is
%   EXACT (slope 1/0), so this is sound; it is also TIGHTER than the no-prefilter path (which widens every
%   bound by derr, so a marginally-stable neuron whose WIDENED CROWN bound straddles 0 gets needlessly
%   relaxed -- the prefilter keeps it exact via the un-widened IBP). NOT bit-identical: a wrong relaxation
%   here = a spuriously-high margin = a wrong unsat = -150, so this asserts the prefilter margin stays a
%   SOUND lower bound on DEEP nets where it genuinely differs from the no-prefilter path.
%     P1 SOUND-vs-ORACLE: prefilter (single) margin <= FP64 oracle margin (the existing sound guarantee).
%     P2 SOUND-vs-MC:     prefilter margin <= the true min of C*f over the box (Monte-Carlo).
%     P3 DIFFERS:         on these deep nets the prefilter genuinely differs from no-prefilter (>=1 of the
%                         set), so P1/P2 are exercising the tightening path, not the trivial identical case.
%   Run: runtests('test_gpu_bab_ibp_prefilter').
    tests = functiontests(localfunctions);
end

function setup(tc)
    tc.TestData.oldSound = getenv('NNV_SOUND_FP32_TIGHT');
    tc.TestData.oldPre   = getenv('NNV_CROWN_IBP_PREFILTER');
    setenv('NNV_SOUND_FP32_TIGHT', '1');
end

function teardown(tc)
    setenv('NNV_SOUND_FP32_TIGHT', tc.TestData.oldSound);
    setenv('NNV_CROWN_IBP_PREFILTER', tc.TestData.oldPre);
end

function test_prefilter_sound_on_deep_nets(tc)
    nin = 12; nout = 6; w = 400;
    C = [1 -1 0 0 0 0; 0 0 1 -1 0 0; 0 0 0 0 1 -1];
    nDiffer = 0; worstOracle = -inf; worstMC = -inf;
    for seed = [1 3 5 9 14 22 31 47]
        rng(seed, 'twister');
        depth = 3 + mod(seed, 4);                          % 3..6 relu layers -> non-trivial derr widening
        ops = {}; pin = nin;
        for L = 1:depth
            ops{end+1} = struct('type','affine','W',randn(w,pin)/sqrt(pin),'b',0.12+0.3*randn(w,1),'src',2*(L-1),'nOut',w); %#ok<AGROW>
            ops{end+1} = struct('type','relu','src',2*L-1); %#ok<AGROW>
            pin = w;
        end
        ops{end+1} = struct('type','affine','W',randn(nout,w)/sqrt(w),'b',randn(nout,1),'src',2*depth,'nOut',nout); %#ok<AGROW>
        lb = -0.1*ones(nin,1); ub = 0.1*ones(nin,1);
        setenv('NNV_CROWN_IBP_PREFILTER', '');  off = double(gpu_bab_crown_tight(ops, lb, ub, C, 'single', {})); off = off(:);
        setenv('NNV_CROWN_IBP_PREFILTER', '1'); on  = double(gpu_bab_crown_tight(ops, lb, ub, C, 'single', {})); on  = on(:);
        setenv('NNV_CROWN_IBP_PREFILTER', '');
        H = double(gpu_bab_crown_tight(ops, lb, ub, C, 'double', {})); H = H(:);
        verifyTrue(tc, all(isfinite(on)), sprintf('seed %d: finite', seed));
        % P1: prefilter single margin <= FP64 oracle (sound). tol 1e-5 >> double roundoff, << a real unsound gap.
        verifyLessThanOrEqual(tc, on, H + 1e-5, sprintf('P1 seed %d: prefilter margin must be <= FP64 oracle', seed));
        worstOracle = max(worstOracle, max(on - H));
        if max(abs(on - off)) > 1e-6
            nDiffer = nDiffer + 1;
            N = 40000; X = lb + (ub - lb).*rand(nin, N); cur = X;
            for L = 1:depth, cur = max(ops{2*L-1}.W*cur + ops{2*L-1}.b, 0); end
            mcmin = min(C*(ops{end}.W*cur + ops{end}.b), [], 2);
            verifyLessThanOrEqual(tc, on, mcmin + 1e-4, sprintf('P2 seed %d: prefilter margin must be <= MC true min', seed));
            worstMC = max(worstMC, max(on - mcmin));
        end
    end
    % P3: the tightening path must actually be exercised (else P1/P2 only tested the trivial identical case).
    verifyGreaterThanOrEqual(tc, nDiffer, 1, 'P3: at least one deep net must differ prefilter-vs-baseline (else not exercising the flip)');
    fprintf('  [ibp-prefilter] differing nets=%d  worst(pre-oracle)=%.2e  worst(pre-MC)=%.2e\n', nDiffer, worstOracle, worstMC);
end
