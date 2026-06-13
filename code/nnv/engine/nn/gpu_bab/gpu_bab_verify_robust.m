function [status, margins] = gpu_bab_verify_robust(ops, x_lb, x_ub, trueLabel, nClasses, precision)
% GPU_BAB_VERIFY_ROBUST  Sound CROWN robustness check (single pass, no branching).
%
%   [status, margins] = GPU_BAB_VERIFY_ROBUST(ops, x_lb, x_ub, trueLabel, nClasses, precision)
%   builds the robustness spec C (rows e_trueLabel - e_j for each other class j) and
%   bounds it with gpu_bab_crown_spec. Returns:
%     status  : 'robust'  -- PROVEN: every input in the box keeps argmax = trueLabel
%               'unknown' -- CROWN bound too loose to prove it (branch-and-bound would
%                            refine; this single pass never returns 'unsafe').
%     margins : (nClasses-1)-by-B lower bounds on out_true - out_j.
%
%   Sound-or-unknown by construction: 'robust' is only returned when a guaranteed
%   lower bound proves every class margin positive; otherwise 'unknown' (never a wrong
%   verdict). Batches over B sub-boxes (margins is per-box). See gpu_bab_bab for the
%   branch-and-bound that turns 'unknown' boxes into 'robust' or a counterexample.

    if nargin < 6 || isempty(precision)
        precision = 'single';
    end
    C = -eye(nClasses, precision);
    C(:, trueLabel) = C(:, trueLabel) + 1;
    C(trueLabel, :) = [];                      % drop the all-zero true-vs-true row
    margins = gpu_bab_crown_spec(ops, x_lb, x_ub, C, precision);
    if all(margins(:) > 0)
        status = 'robust';
    else
        status = 'unknown';
    end
end
