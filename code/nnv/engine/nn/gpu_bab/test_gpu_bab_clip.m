function tests = test_gpu_bab_clip
% TEST_GPU_BAB_CLIP  Soundness of Clip-and-Verify Relaxed Clipping (gpu_bab_clip):
%   the clipped box must CONTAIN the constraint-feasible subset of the original box (sound
%   over-approx), must TIGHTEN when a constraint cuts a coordinate, and must report EMPTY only
%   when the constraint system is genuinely infeasible over the box (no false "empty" -> no
%   wrong UNSAT). Run: runtests('test_gpu_bab_clip').
    tests = functiontests(localfunctions);
end

function test_clip_tighten_and_contain(tc)
    rng(3, 'twister');
    n = 4; lb = -ones(n,1); ub = ones(n,1);
    a = [1;0;0;0]; c = -0.3;                 % constraint x1 <= 0.3
    [lb2, ub2, empty] = gpu_bab_clip(lb, ub, a.', c);
    verifyFalse(tc, empty, 'feasible constraint wrongly reported empty');
    verifyEqual(tc, ub2(1), 0.3, 'AbsTol', 1e-9);    % coord 1 tightened to 0.3
    verifyEqual(tc, lb2(1), -1, 'AbsTol', 1e-9);     % lower unchanged
    % containment: every box point satisfying the constraint lies in [lb2,ub2]
    N = 30000; X = lb + (ub - lb) .* rand(n, N);
    Xf = X(:, (a.' * X + c) <= 0);
    verifyTrue(tc, all(Xf >= lb2 - 1e-9, 'all') && all(Xf <= ub2 + 1e-9, 'all'), ...
        'clipped box does not contain the feasible set (UNSOUND)');
end

function test_clip_multi_constraint_contain(tc)
    rng(4, 'twister');
    n = 3; lb = [-1;-1;-1]; ub = [1;1;1];
    A = [1 1 0; 0 1 1]; c = [-0.4; -0.6];    % x1+x2<=0.4 AND x2+x3<=0.6
    [lb2, ub2, empty] = gpu_bab_clip(lb, ub, A, c);
    verifyFalse(tc, empty);
    N = 40000; X = lb + (ub - lb) .* rand(n, N);
    feas = all(A * X + c <= 0, 1);
    Xf = X(:, feas);
    if ~isempty(Xf)
        verifyTrue(tc, all(Xf >= lb2 - 1e-9, 'all') && all(Xf <= ub2 + 1e-9, 'all'), ...
            'multi-constraint clipped box does not contain the feasible set (UNSOUND)');
    end
end

function test_clip_detects_empty(tc)
    n = 4; lb = -ones(n,1); ub = ones(n,1);
    % single infeasible constraint: x1 <= -2 over x1 in [-1,1]
    [~,~,e1] = gpu_bab_clip(lb, ub, [1 0 0 0], -2);   % x1 + (-2) ... a'x+c = x1-2 <= 0 always true? no
    % a'x + c <= 0 with a=[1..],c=-2 is x1 <= 2 -> ALWAYS feasible; use c=+2 -> x1 <= -2 infeasible
    [~,~,e1b] = gpu_bab_clip(lb, ub, [1 0 0 0], 2);
    verifyTrue(tc, e1b, 'x1<=-2 over [-1,1] should be empty');
    verifyFalse(tc, e1, 'x1<=2 over [-1,1] should be feasible');
    % contradictory pair: x1>=0.6 (-x1+0.6<=0) AND x1<=-0.6 (x1+0.6<=0)
    A = [-1 0 0 0; 1 0 0 0]; c = [0.6; 0.6];
    [~,~,e2] = gpu_bab_clip(lb, ub, A, c);
    verifyTrue(tc, e2, 'contradictory constraints should be empty');
end
