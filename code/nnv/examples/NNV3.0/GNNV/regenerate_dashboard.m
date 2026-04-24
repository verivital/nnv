% regenerate_dashboard - Regenerate dashboard and table from saved results
% Useful for updating figures without re-running verification experiments

% Load saved results
load('results/gnn_results.mat', 'results');
load('models/gine_ieee24.mat', 'src', 'dst');
model_data = struct('src', src, 'dst', dst);
figuresDir = 'figures';

% Regenerate dashboard and table
generate_cav26_dashboard(results, model_data, figuresDir);
generate_latex_table(results, figuresDir);

fprintf('Dashboard regenerated from saved results!\n');
