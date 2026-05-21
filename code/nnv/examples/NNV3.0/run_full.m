% run_full.m -- ~3-5 h full reproduction of the NNV3.0 artifact in the
% CURRENT MATLAB session.
%
% Same as run_smoke.m but ToolComparison runs in 'full' mode and
% renders the ATVA 2026 paper's Tables 5, 6, and 7 to
% ToolComparison/tables/out/table_main.{tex,txt}.
%
% Usage (paste in the MATLAB Command Window after setup_online):
%   run('/home/matlab/nnv/code/nnv/examples/NNV3.0/run_full.m')

addpath('/home/matlab/nnv/code/nnv/examples/NNV3.0/utils');
run_experiments('full');
