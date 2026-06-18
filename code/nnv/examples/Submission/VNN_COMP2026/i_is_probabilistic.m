function tf = i_is_probabilistic(reachMethod)
% I_IS_PROBABILISTIC  TRUE iff a reach method's 'holds'/unsat is NOT a sound proof.
%   cp-star (Prob_reach, conformal prediction) is the only probabilistic reach method in the
%   VNN-COMP harness: its verdict is a statistical guarantee, not a sound over-approximation, so a
%   cp-star 'unsat' can be wrong (-150 in VNN-COMP scoring). Every sound method (approx-star,
%   relax-star-area, exact-star, approx-zono, abs-dom) returns FALSE.
%
%   Centralized in ONE place (its own file, so it is unit-testable) and used by BOTH the cp-star
%   quarantine (i_finalize_reach_options) and the cp-star verdict log (i_log_cpstar), so the
%   definition of "probabilistic" can never drift between the gate and the label.
%
%   See also: i_finalize_reach_options, run_vnncomp_instance,
%             tests/nn/vnncomp/test_run_vnncomp_routing.m
    tf = strcmp(reachMethod, 'cp-star');
end
