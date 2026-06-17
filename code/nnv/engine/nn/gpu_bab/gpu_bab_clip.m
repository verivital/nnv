function [lb2, ub2, empty] = gpu_bab_clip(lb, ub, A, c)
% GPU_BAB_CLIP  Clip-and-Verify Relaxed Clipping (Zhou et al., NeurIPS 2025, Thm 3.2).
%   [lb2, ub2, empty] = GPU_BAB_CLIP(lb, ub, A, c) tightens the axis-aligned box [lb,ub] against
%   the linear constraint system A*x + c <= 0 (A is m-by-n, c is m-by-1). It returns the tightest
%   box that still contains the feasible set F = { x in [lb,ub] : A*x + c <= 0 }, and a flag
%   `empty` that is true when F is provably EMPTY (the constraint system is infeasible over the
%   box). One O(n) closed-form pass per constraint (parallel version: every constraint clips from
%   the ORIGINAL box center, then the elementwise-tightest bound is kept -- order-independent and
%   a sound over-approximation).
%
%   USE FOR JOINT POLYTOPE AVOIDANCE (the gap single-halfspace separation misses): to prove the
%   network output avoids an unsafe disjunct {G*y <= g}, take the CROWN input-space LOWER plane of
%   each row, G_i*y >= Ain_i*x + din_i (from gpu_bab_crown_tight), so the unsafe set
%   {x : G_i*y(x) <= g_i for ALL i} is a SUBSET of {x : Ain_i*x + (din_i - g_i) <= 0 for all i}.
%   Call gpu_bab_clip(lb, ub, Ain, din - g); if `empty`, NO input reaches the unsafe polytope ->
%   the disjunct is avoided JOINTLY -> sound UNSAT. (Each row alone need not separate.)
%
%   SOUNDNESS (Thm 3.2 + its proof B.4): for a single constraint a'x + c <= 0, the per-coordinate
%   clip x_i^clip = ( -sum_{j!=i}[a_j*xhat_j - |a_j|*eps_j] - c ) / a_i (xhat = box center, eps =
%   half-width) is the exact max/min of x_i over F, so F subset [lb2,ub2] subset [lb,ub]. The
%   infeasibility test min_{box}(a'x + c) = a'xhat - sum|a|eps + c > 0 means the constraint cannot
%   hold anywhere in the box. For multiple constraints, F subset (intersection of per-constraint
%   over-approx boxes) subset [lb2,ub2]; empty => F empty. No false "empty" -> no wrong UNSAT.

    lb = double(lb(:)); ub = double(ub(:)); n = numel(lb);
    A = double(A); c = double(c(:)); m = size(A, 1);
    lb2 = lb; ub2 = ub; empty = false;
    if m == 0, return; end
    if size(A, 2) ~= n
        error('gpu_bab_clip:dim', 'A has %d cols but the box is %d-dim.', size(A,2), n);
    end
    xhat = (lb + ub) / 2;
    eps  = (ub - lb) / 2;                          % >= 0
    tol  = 1e-12;
    for k = 1:m
        a = A(k, :).';                             % n x 1
        ck = c(k);
        boxMin = a.' * xhat - sum(abs(a) .* eps);  % min over the box of a'x
        if boxMin + ck > tol                       % constraint a'x + c <= 0 infeasible over the box
            empty = true; return;
        end
        nz = find(a ~= 0);
        for q = 1:numel(nz)
            i = nz(q); ai = a(i);
            % worst (most-negative) value of sum_{j != i} a_j x_j over the box:
            crossMin = boxMin - (ai * xhat(i) - abs(ai) * eps(i));
            xclip = (-crossMin - ck) / ai;
            if ai > 0
                ub2(i) = min(ub2(i), xclip);
            else
                lb2(i) = max(lb2(i), xclip);
            end
        end
    end
    if any(lb2 > ub2 + tol)                        % some coordinate clipped to empty
        empty = true;
    end
end
