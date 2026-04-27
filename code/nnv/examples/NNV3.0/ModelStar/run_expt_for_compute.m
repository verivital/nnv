% ModelStar — perturb weights of single layers in an MNIST MLP and verify
% robustness under each perturbation magnitude. Writes per-fraction
% results into a freshly built experiment configuration (results/MNIST_MLP.mat).
%
% Run: matlab -batch "cd code/nnv/examples/NNV3.0/ModelStar; run_expt_for_compute"
%
% Caller may pre-set `n_layers_to_run_for_from_yaml_file` (e.g., =3 to
% reproduce the paper's fc_4/fc_5/fc_6 sweep, or =1 for a quick sanity
% run on fc_6 only); the default applies otherwise.

% Run from this script's directory so relative paths into results/ work.
ms_scriptDir = fileparts(mfilename('fullpath'));
cd(ms_scriptDir);

if ~exist('results', 'dir'); mkdir('results'); end
results_file = fullfile('results', 'MNIST_MLP');

% Build the empty experiment template programmatically (replaces the
% legacy YAML template — no external yaml package required).
data = build_template(); %#ok<NASGU>  consumed by save below
builtin('save', [results_file '.mat'], 'data');

do_not_clear = 1;
if ~exist('n_layers_to_run_for_from_yaml_file', 'var')
    % Paper sweep: fc_6, fc_5, fc_4 (the last three layers).
    n_layers_to_run_for_from_yaml_file = 3;
end
expt_path = results_file;
conv_expt_any_layer
