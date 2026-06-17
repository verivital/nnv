function tests = test_gpu_bab_halfspace_genbab
% TEST_GPU_BAB_HALFSPACE_GENBAB  Soundness of the GenBaB bilinear halfspace verifier
%   (gpu_bab_halfspace_genbab). Builds a ReLU + bilinear-product net (the lsnc Lyapunov shape:
%   input -> affine -> ReLU -> affine, and input -> affine, multiplied, -> affine spec output) and
%   checks the two properties that matter for a verdict-emitting verifier:
%     (1) SOUND-OR-UNKNOWN: when the unsafe region is genuinely reachable (Monte-Carlo finds an
%         input whose output lands inside it), the verifier must NEVER answer 'robust'. A wrong
%         'robust' is a -150 VNN-COMP penalty, so this is the load-bearing test.
%     (2) CERTIFIES a clearly-separated spec: when the unsafe region sits far outside the reachable
%         set, the verifier answers 'robust' (the product + ReLU branching + Clip-and-Verify path
%         actually closes a safe instance, not just bails to 'unknown').
%   Run: runtests('test_gpu_bab_halfspace_genbab').
    tests = functiontests(localfunctions);
end

function net = i_net()
    rng(7, 'twister');
    W1 = randn(5, 3); b1 = randn(5, 1);     % op1: affine 3->5
    W2 = randn(3, 5); b2 = randn(3, 1);     % op3: affine 5->3  (product input 1, after ReLU)
    W3 = randn(3, 3); b3 = randn(3, 1);     % op4: affine 3->3  (product input 2, from raw input)
    W4 = randn(1, 3); b4 = randn(1, 1);     % op6: affine 3->1  (scalar spec output)
    net.ops = { struct('type','affine','W',W1,'b',b1,'src',0), ...
                struct('type','relu','src',1), ...
                struct('type','affine','W',W2,'b',b2,'src',2), ...
                struct('type','affine','W',W3,'b',b3,'src',0), ...
                struct('type','product','inputs',[3 4],'sizes',[3 3],'src',3), ...
                struct('type','affine','W',W4,'b',b4,'src',5) };
    net.lb = [-0.2; -0.2; -0.2];
    net.ub = [ 0.2;  0.2;  0.2];
end

function [sMin, sMax] = i_mc_range(ops, lb, ub)
    rng(3, 'twister');
    N = 80000;
    X = lb + (ub - lb) .* rand(numel(lb), N);
    S = i_eval(ops, X);            % 1 x N scalar spec output
    sMin = min(S); sMax = max(S);
end

function test_sound_when_unsafe(tc)
    % Unsafe region {s <= sMid} is reachable (~half the box) -> verifier must NOT say 'robust'.
    net = i_net();
    [sMin, sMax] = i_mc_range(net.ops, net.lb, net.ub);
    sMid = (sMin + sMax) / 2;
    Gd = {1}; gd = {sMid};                      % unsafe disjunct: 1*s <= sMid
    % sanity: the unsafe region is genuinely reachable
    rng(5, 'twister'); X = net.lb + (net.ub - net.lb) .* rand(3, 40000);
    S = i_eval(net.ops, X);
    verifyTrue(tc, any(S <= sMid), 'test setup: unsafe region should be reachable');
    opts = struct('maxNodes', 6000, 'timeCap', 30, 'precision', 'double');
    [verdict, ~] = gpu_bab_halfspace_genbab(net.ops, net.lb, net.ub, Gd, gd, opts);
    verifyFalse(tc, strcmp(verdict, 'robust'), ...
        'UNSOUND: returned robust for a reachable unsafe region (-150 risk)');
end

function test_certifies_when_separated(tc)
    % Unsafe region {s <= sMin - 20} sits far below the reachable set -> certifiable 'robust'.
    net = i_net();
    [sMin, ~] = i_mc_range(net.ops, net.lb, net.ub);
    Gd = {1}; gd = {sMin - 20};
    opts = struct('maxNodes', 50000, 'timeCap', 60, 'precision', 'double');
    [verdict, info] = gpu_bab_halfspace_genbab(net.ops, net.lb, net.ub, Gd, gd, opts);
    verifyEqual(tc, verdict, 'robust', ...
        sprintf('expected robust on a well-separated spec (reason: %s)', info.reason));
end

function test_disjunction_sound(tc)
    % Two-disjunct unsafe set, both reachable -> must NOT say 'robust'.
    net = i_net();
    [sMin, sMax] = i_mc_range(net.ops, net.lb, net.ub);
    Gd = {1; -1}; gd = {sMax; -(sMin)};         % {s<=sMax} OR {-s<=-sMin} i.e. {s>=sMin}: covers all
    opts = struct('maxNodes', 4000, 'timeCap', 20, 'precision', 'double');
    [verdict, ~] = gpu_bab_halfspace_genbab(net.ops, net.lb, net.ub, Gd, gd, opts);
    verifyFalse(tc, strcmp(verdict, 'robust'), ...
        'UNSOUND: returned robust though the union of unsafe disjuncts covers the reachable set');
end

function test_fuzz_never_false_robust(tc)
    % Adversarial fuzz for the off-region McCormick concern (soundness review finding 2): across many
    % random bilinear nets and VIOLATED specs (unsafe region IS reachable, MC-confirmed, deliberately
    % SMALL so a violation tends to sit inside ONE product-split sibling region), the verifier must
    % NEVER return 'robust'. The product-input override makes the McCormick plane valid only over a
    % child's value-region; if certify could lean on that plane OFF that region, some trial here would
    % slip through as a false 'robust'. Empirically confirms the region-validity soundness argument.
    nTrials = 25; nViolated = 0;
    for trial = 1:nTrials
        rng(1000 + trial, 'twister');
        n = 3;
        W1 = randn(5, n); b1 = randn(5, 1);
        W2 = randn(3, 5); b2 = randn(3, 1);
        W3 = randn(3, n); b3 = randn(3, 1);
        W4 = randn(1, 3); b4 = randn(1, 1);
        ops = { struct('type','affine','W',W1,'b',b1,'src',0), ...
                struct('type','relu','src',1), ...
                struct('type','affine','W',W2,'b',b2,'src',2), ...
                struct('type','affine','W',W3,'b',b3,'src',0), ...
                struct('type','product','inputs',[3 4],'sizes',[3 3],'src',3), ...
                struct('type','affine','W',W4,'b',b4,'src',5) };
        lb = -0.2 * ones(n, 1); ub = 0.2 * ones(n, 1);
        N = 30000; X = lb + (ub - lb) .* rand(n, N); S = i_eval(ops, X);
        sMin = min(S); sMax = max(S);
        g = sMin + 0.08 * (sMax - sMin);          % unsafe = {s<=g}: reachable but small (~bottom 8%)
        if ~any(S <= g), continue; end            % keep only genuinely violated trials
        nViolated = nViolated + 1;
        opts = struct('maxNodes', 2500, 'timeCap', 15, 'precision', 'double');
        [verdict, ~] = gpu_bab_halfspace_genbab(ops, lb, ub, {1}, {g}, opts);
        verifyFalse(tc, strcmp(verdict, 'robust'), ...
            sprintf('UNSOUND on trial %d: robust returned for a reachable (violated) unsafe region', trial));
    end
    verifyGreaterThan(tc, nViolated, 15, 'fuzz setup: too few genuinely-violated trials to be meaningful');
end

function Y = i_eval(ops, X)
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
