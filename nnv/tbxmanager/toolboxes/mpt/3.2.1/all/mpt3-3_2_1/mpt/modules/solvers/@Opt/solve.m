function [ret,opt] = solve(opt, varargin)
%
% Solve the problem
%

global MPTOPTIONS

if isempty(MPTOPTIONS)
    MPTOPTIONS = mptopt;
end

% deal with arrays
ret=cell(size(opt));
if numel(opt)>1
    for i=1:numel(opt)
        ret{i} = opt(i).solve;
    end
    return
end

% Determine problem type
if opt.isParametric && nargin==1
    % for parametric solvers use "mpt_solvemp" which depends
    % on OPT and POLYHEDRON classes
    ret = mpt_solvemp(opt);
else
    % non-parametric solvers can be called directly.
	%
	% alternatively, we can solve a parametric problem for a particular
	% value of the parameter
    ret = mpt_solve(opt, varargin{:});
end

end
