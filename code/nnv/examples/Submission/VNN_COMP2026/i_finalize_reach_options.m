function reachOptionsList = i_finalize_reach_options(reachOptionsList, category, nnvnetValid)
% I_FINALIZE_REACH_OPTIONS  Pure method-selection post-processing for the VNN-COMP harness.
%   Extracted from load_vnncomp_network (run_vnncomp_instance.m) into its own file so the routing
%   -- i.e. WHICH reach methods actually run per category, and whether a PROBABILISTIC method can
%   survive -- is guarded by unit tests instead of being implicit in a ~1500-line function. It is
%   pure (no MATLAB import, no GPU, no reach), so a test can call it directly with a synthetic
%   reachOptionsList and assert the resulting method order.
%
%   Applies, in order, exactly what the inline code did, plus a new opt-in quarantine:
%     (1) Phase-1.5 sound prepend -- if the list LEADS with cp-star (probabilistic) AND the net is
%         valid, prepend {approx-star, relax-star-area 0.5} so the SOUND methods are tried first;
%         cp-star stays in the list as the last-resort fallback.
%     (2) slow_cats {cifar100, cora, safenlp, sat_relu, tinyimagenet, vggnet} -- DROP exact-star and
%         lead with the fast, looser, but STILL-SOUND over-approximations {approx-zono, abs-dom} so
%         a compute-bound category produces a sound verdict instead of timing out.
%     (3) cp-star QUARANTINE -- gated by env NNV_QUARANTINE_CPSTAR (DEFAULT OFF). When set, STRIP
%         cp-star entirely so a probabilistic 'unsat' can never be emitted as a sound verdict in ANY
%         reach branch (single-spec OR either parfor branch). Strictly sound: removing cp-star can
%         only turn a would-be unsat into 'unknown', never produce a wrong verdict. DEFAULT OFF =>
%         the harness behaves exactly as before; set the env to measure / enforce sound-or-unknown.
%
%   Inputs:
%     reachOptionsList : 1xN cell of reachOptions structs (each with a .reachMethod field)
%     category         : the VNN-COMP category string (substring-matched, e.g. 'cifar100_2024')
%     nnvnetValid      : logical, is_nnvnet_valid(nnvnet) -- both prepends are gated on it, matching
%                        the original code (a sentinel net "" disables the prepends).
%
%   See also: i_is_probabilistic, is_nnvnet_valid, run_vnncomp_instance,
%             tests/nn/vnncomp/test_run_vnncomp_routing.m
    if isempty(reachOptionsList)
        return;
    end

    % (1) sound prepend ahead of a cp-star-led list
    if isfield(reachOptionsList{1}, 'reachMethod') ...
            && strcmp(reachOptionsList{1}.reachMethod, 'cp-star') && nnvnetValid
        o1 = struct(); o1.reachMethod = 'approx-star';
        o2 = struct(); o2.reachMethod = 'relax-star-area'; o2.relaxFactor = 0.5;
        reachOptionsList = [{o1, o2}, reachOptionsList];
    end

    % (2) slow_cats: drop exact-star, lead with zonotope + abstract-domain
    slow_cats = ["cifar100","cora","safenlp","sat_relu","tinyimagenet","vggnet"];
    if any(contains(category, slow_cats)) && nnvnetValid
        kept = reachOptionsList(~cellfun(@(o) strcmp(o.reachMethod, 'exact-star'), reachOptionsList));
        zo = struct(); zo.reachMethod = 'approx-zono';
        ad = struct(); ad.reachMethod = 'abs-dom';
        reachOptionsList = [{zo, ad}, kept];
    end

    % (3) cp-star quarantine (env-gated, DEFAULT OFF -> identity / behavior unchanged)
    if ~isempty(getenv('NNV_QUARANTINE_CPSTAR'))
        reachOptionsList = reachOptionsList(~cellfun(@(o) i_is_probabilistic(o.reachMethod), reachOptionsList));
    end
end
