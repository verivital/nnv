function sol = fast_chebyCenter(H, He, Rmax)
% Computes the center and the radius of a chebyball
%
%    sol = fast_chebyCenter(H, He, Rmax)
%
% For internal purposes only! Use Polyhedron/chebyCenter() instead

global MPTOPTIONS

narginchk(3, 3);
A = H(:, 1:end-1);
S.A = [A sqrt(sum(A.*A,2))];
S.b = H(:, end);
S.Ae = [He(:, 1:end-1) zeros(size(He,1),1)];
S.be = He(:, end);
% lower bounds on [xc, r]
S.A = [S.A; zeros(1,size(A, 2)), -1];
S.b = [S.b; 0];
% upper bound on the chebyradius
if ~isinf(Rmax)
    S.A = [S.A; zeros(1, size(A, 2)), 1];
    S.b = [S.b; Rmax];
end
% the last value is -1 because it is maximization
S.f = [zeros(size(A, 2), 1); -1];
% solve the cheby problem
S.lb = []; S.ub = []; S.quicklp = true;
ret = mpt_solve(S);
% set output variables
sol.exitflag = ret.exitflag;
if ret.exitflag == MPTOPTIONS.OK
	if -ret.obj>MPTOPTIONS.zero_tol
		sol.x = ret.xopt(1:end-1);
		sol.r = ret.xopt(end);
	else
		sol.x = ret.xopt(1:end-1);
		sol.r = 0;
	end
else
	% infeasible
	sol.x = [];
	sol.r = -Inf;
end

end
