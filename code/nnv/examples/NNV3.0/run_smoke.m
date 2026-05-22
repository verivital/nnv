% run_smoke.m -- ~20-25 min smoke test of the NNV3.0 artifact, run in
% the CURRENT MATLAB session.
%
% MATLAB equivalent of `bash run_all.sh` for the nnv3.0-online flow:
% reuses this session's browser sign-in licence rather than spawning
% per-experiment matlab -batch processes (which would each need their
% own licence checkout, and online sign-in only works for -browser).
%
% Usage (paste in the MATLAB Command Window after setup_online):
%   run('/home/matlab/nnv/code/nnv/examples/NNV3.0/run_smoke.m')
%
% ToolComparison runs in 'smoke' mode (~12 min, NNV-only sanity pass).
% Use run_full.m for the ~3-5 h Tables-5/6/7 reproduction.

addpath('/home/matlab/nnv/code/nnv/examples/NNV3.0/utils');
run_experiments('smoke');
