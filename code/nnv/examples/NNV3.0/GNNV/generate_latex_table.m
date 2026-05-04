function generate_latex_table(results, figuresDir)
% generate_latex_table - Generate LaTeX table with verification percentages and timing
%
% Creates a publication-ready LaTeX table showing verification percentages
% and average verification time across all models and epsilon values.
% Uses booktabs style for professional appearance.
%
% Usage:
%   generate_latex_table(results, figuresDir)
%
% Inputs:
%   results    - Results struct from run_gnn_experiments
%   figuresDir - Output directory for the .tex file
%
% Output:
%   figures/results_table.tex
%
% Author: Anne Tumlin
% Date: 01/21/2026

models = results.config.models;
epsilons = results.config.epsilons;
num_scenarios = results.config.num_scenarios;

% Total nodes for percentage calculation
total_nodes = results.data{1,1,1}.verified + results.data{1,1,1}.unknown + results.data{1,1,1}.violated;
total_possible = total_nodes * num_scenarios;  % 13 * 10 = 130

% Calculate total percentages and average times
pct_total = zeros(length(models), length(epsilons));
avg_time = zeros(length(models), length(epsilons));

for m = 1:length(models)
    for e = 1:length(epsilons)
        total_verified = 0;
        total_time = 0;
        for s = 1:num_scenarios
            total_verified = total_verified + results.data{m, e, s}.verified;
            total_time = total_time + results.data{m, e, s}.time;
        end
        pct_total(m, e) = (total_verified / total_possible) * 100;
        avg_time(m, e) = total_time / num_scenarios;
    end
end

% Build LaTeX table with booktabs style
fid = fopen(fullfile(figuresDir, 'results_table.tex'), 'w');

fprintf(fid, '\\begin{table}[h]\n');
fprintf(fid, '\\centering\n');
fprintf(fid, '\\caption{Verification results for GNN-based power flow prediction on the IEEE 24-bus system. ');
fprintf(fid, 'Percentage of voltage nodes verified safe (out of 130 total across 10 scenarios) ');
fprintf(fid, 'and average time to verify all 13 voltage nodes per graph.}\n');
fprintf(fid, '\\label{tab:gnn-results}\n');
fprintf(fid, '\\begin{tabular}{l @{\\hspace{0.8em}} cc @{\\hspace{0.8em}} cc @{\\hspace{0.8em}} cc}\n');
fprintf(fid, '\\toprule\n');

% Header row 1: epsilon values spanning two columns each
fprintf(fid, '& \\multicolumn{2}{c}{$\\epsilon=%.3f$} & \\multicolumn{2}{c}{$\\epsilon=%.3f$} & \\multicolumn{2}{c}{$\\epsilon=%.2f$} \\\\\n', ...
    epsilons(1), epsilons(2), epsilons(3));

% Add cmidrule for visual separation under epsilon headers
fprintf(fid, '\\cmidrule(lr){2-3} \\cmidrule(lr){4-5} \\cmidrule(lr){6-7}\n');

% Header row 2: % and Time labels
fprintf(fid, 'Model & Verified (\\%%) & Time (s) & Verified (\\%%) & Time (s) & Verified (\\%%) & Time (s) \\\\\n');
fprintf(fid, '\\midrule\n');

% Data rows
for m = 1:length(models)
    fprintf(fid, '%s', models{m});
    for e = 1:length(epsilons)
        fprintf(fid, ' & %.1f & %.2f', pct_total(m, e), avg_time(m, e));
    end
    fprintf(fid, ' \\\\\n');
end

fprintf(fid, '\\bottomrule\n');
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\end{table}\n');

fclose(fid);

fprintf('LaTeX table saved to: %s\n', fullfile(figuresDir, 'results_table.tex'));
end
