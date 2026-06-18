%% test_run_vnncomp_routing
% Routing / method-selection regression tests for the VNN-COMP harness (run_vnncomp_instance).
%
% WHY THIS EXISTS: a missing test let a real soundness fragility slip -- the harness can fall
% through to cp-star (Prob_reach = conformal prediction = PROBABILISTIC) and emit its 'unsat' as
% if it were a sound proof. These tests pin down WHICH reach methods actually run per category
% (the thing we could not previously tell from a result) and that the NNV_QUARANTINE_CPSTAR gate
% removes every probabilistic method. They exercise the PURE, extracted routing helpers
% (i_is_probabilistic, i_finalize_reach_options) directly -- no ONNX import, no GPU, no reach --
% so they are fast, deterministic, and run anywhere.
%
% Run: runtests('test_run_vnncomp_routing')
%
% Each section sets its OWN env + path so the sections are order-independent (the script-test
% framework runs %% sections with independent workspaces; env + path are process-global, so we
% set them explicitly per section rather than relying on a prior section).

%% i_is_probabilistic identifies ONLY cp-star
addpath(fileparts(mfilename('fullpath')));
assert(i_is_probabilistic('cp-star') == true,  'cp-star must be probabilistic');
assert(i_is_probabilistic('approx-star') == false);
assert(i_is_probabilistic('relax-star-area') == false);
assert(i_is_probabilistic('exact-star') == false);
assert(i_is_probabilistic('approx-zono') == false);
assert(i_is_probabilistic('abs-dom') == false);

%% cifar100 effective ladder = 4 SOUND rungs THEN cp-star LAST (default OFF, valid net)
addpath(fileparts(mfilename('fullpath')));
setenv('NNV_QUARANTINE_CPSTAR','');                       % default behavior
L = i_finalize_reach_options(local_mkopts({'cp-star'}), 'cifar100_2024', true);
names = cellfun(@(o) o.reachMethod, L, 'uni', 0);
assert(isequal(names, {'approx-zono','abs-dom','approx-star','relax-star-area','cp-star'}), ...
    ['cifar100 ladder wrong: ' strjoin(names, ',')]);
assert(i_is_probabilistic(names{end}) == true, 'cp-star must be the LAST (last-resort) rung');
assert(all(cellfun(@(n) ~i_is_probabilistic(n), names(1:end-1))), 'every leading rung must be SOUND');

%% QUARANTINE strips cp-star -> NO probabilistic method survives for cifar100 (sound rungs kept)
addpath(fileparts(mfilename('fullpath')));
setenv('NNV_QUARANTINE_CPSTAR','1');
L = i_finalize_reach_options(local_mkopts({'cp-star'}), 'cifar100_2024', true);
names = cellfun(@(o) o.reachMethod, L, 'uni', 0);
assert(~any(cellfun(@i_is_probabilistic, names)), 'quarantine must remove cp-star');
assert(isequal(names, {'approx-zono','abs-dom','approx-star','relax-star-area'}), 'sound rungs must remain');
setenv('NNV_QUARANTINE_CPSTAR','');

%% ml4acopf-style cp-star-ONLY on an INVALID net is the -150 path; quarantine -> empty (sound)
addpath(fileparts(mfilename('fullpath')));
setenv('NNV_QUARANTINE_CPSTAR','');
Lraw = i_finalize_reach_options(local_mkopts({'cp-star'}), 'ml4acopf', false);   % invalid net -> no prepend
assert(isequal(cellfun(@(o) o.reachMethod, Lraw, 'uni',0), {'cp-star'}), ...
    'documents the landmine: cp-star ALONE on an invalid net');
setenv('NNV_QUARANTINE_CPSTAR','1');
Lq = i_finalize_reach_options(local_mkopts({'cp-star'}), 'ml4acopf', false);
assert(isempty(Lq), 'quarantine of a cp-star-only list -> empty -> falsify/SAT-or-unknown only (sound)');
setenv('NNV_QUARANTINE_CPSTAR','');

%% cersyve trailing cp-star: quarantine removes ONLY cp-star, keeps the two sound rungs
addpath(fileparts(mfilename('fullpath')));
setenv('NNV_QUARANTINE_CPSTAR','1');
L = i_finalize_reach_options(local_mkopts({'approx-star','relax-star-area','cp-star'}), 'cersyve', true);
names = cellfun(@(o) o.reachMethod, L, 'uni', 0);
assert(isequal(names, {'approx-star','relax-star-area'}), 'cersyve: cp-star dropped, sound rungs intact');
assert(~any(cellfun(@i_is_probabilistic, names)));
setenv('NNV_QUARANTINE_CPSTAR','');

%% a SOUND-only category (acasxu exact-star) is untouched by every step (incl. quarantine)
addpath(fileparts(mfilename('fullpath')));
setenv('NNV_QUARANTINE_CPSTAR','1');                      % even with quarantine on, no cp-star to strip
L = i_finalize_reach_options(local_mkopts({'exact-star'}), 'acasxu', true);
assert(isequal(cellfun(@(o) o.reachMethod, L, 'uni',0), {'exact-star'}), 'acasxu must be unchanged');
setenv('NNV_QUARANTINE_CPSTAR','');

%% slow_cats DROP exact-star and LEAD with zono/abs-dom (timeout guidance), valid net
addpath(fileparts(mfilename('fullpath')));
setenv('NNV_QUARANTINE_CPSTAR','');
L = i_finalize_reach_options(local_mkopts({'approx-star','exact-star'}), 'safenlp', true);
names = cellfun(@(o) o.reachMethod, L, 'uni', 0);
assert(~any(strcmp(names,'exact-star')), 'slow_cats must drop exact-star');
assert(isequal(names(1:2), {'approx-zono','abs-dom'}), 'slow_cats must lead with zono/abs-dom');

%% INVALID net disables BOTH prepends (matches is_nnvnet_valid gating in the original code)
addpath(fileparts(mfilename('fullpath')));
setenv('NNV_QUARANTINE_CPSTAR','');
L = i_finalize_reach_options(local_mkopts({'cp-star'}), 'cifar100_2024', false);
assert(isequal(cellfun(@(o) o.reachMethod, L, 'uni',0), {'cp-star'}), 'invalid net -> no prepend');

%% empty input list is returned unchanged (cgan/nn4sys-mscn falsify-only categories)
addpath(fileparts(mfilename('fullpath')));
setenv('NNV_QUARANTINE_CPSTAR','1');
L = i_finalize_reach_options({}, 'cgan2026', true);
assert(isempty(L), 'empty reachOptionsList stays empty (falsify-only is sound)');
setenv('NNV_QUARANTINE_CPSTAR','');

%% Summary (trailing deterministic section; avoids the script-test line-number quirk)
disp('test_run_vnncomp_routing: all routing/method-selection assertions passed');

% ---- local helpers (file-scope; available to every section above) ----
function L = local_mkopts(methods)
    % Build a reachOptionsList (1xN cell of structs with a .reachMethod field) from a cellstr of
    % method names -- the loop lives INSIDE this function (not in a %% section) per the script-test
    % line-number quirk.
    L = cell(1, numel(methods));
    for k = 1:numel(methods)
        o = struct(); o.reachMethod = methods{k};
        L{k} = o;
    end
end
