function ret = distance(P, S)
% DISTANCE Compute the distance between the given point/polyhedron
% and this polyhedron. 
%
% Solves the optimization problem:
%
%   min ||x - y||_2 s.t. y in P, x in S
%
% Parameters:
%   P - Polyhedron
%   S - Polyhedron or a point given as vector
%
% Returns:
%  ret - Structure containing the elements:
%    dist - Distance between S and this Polyhedron
%    y   - Closest element y in this polyhedron
%    x   - Closest element x in S
%  or [] is either set is empty
%

global MPTOPTIONS
if isempty(MPTOPTIONS)
    MPTOPTIONS = mptopt;
end

assert(isa(P, 'Polyhedron'), 'The first input must be a Polyhedron object.');

% check if S is an array or polyhedron
if isa(S,'Polyhedron')
    if numel(S)>1
       error('Only single polyhedron S allowed.'); 
    end
else
    validate_realvector(S);
end

if numel(P)>1
	% return an array of structures
	ret = P.forEach(@(elem) elem.distance(S));
    return
end

% pre-alocate output
ret = struct('exitflag', [], 'dist', Inf, 'x', [], 'y', []);

%% S is supposedly a point
if ~isa(S,'Polyhedron') 
  % Call the superclass, which handles this case
  ret = P.distance@ConvexSet(S);
  return
end

%% S is Polyhedron
if S.Dim ~= P.Dim
    error('Both polyhedra have to be of the same dimension.');
end

if P.isEmptySet() || S.isEmptySet()
    % distance from an empty set is infinite by convention (issue #111)
    ret.exitflag = MPTOPTIONS.INFEASIBLE;
    ret.dist = Inf;
    return
end

% Get representations of both polyhedra
matP = P.optMat;
matS = S.optMat;

% Build optimization matrices
qp.A  = blkdiag(matP.A, matS.A);
qp.b  = [matP.b;matS.b];
qp.Ae = blkdiag(matP.Ae, matS.Ae);
qp.be = [matP.be; matS.be];
qp.lb = [matP.lb; matS.lb];
qp.ub = [matP.ub; matS.ub];

% Build the cost function (x-y)'(x-y)
nLamP = size(matP.A,2) - P.Dim;
nLamS = size(matS.A,2) - P.Dim;
HP  = diag([ones(P.Dim,1);zeros(nLamP,1)]);
HPS = blkdiag(-eye(P.Dim), zeros(nLamP, nLamS));
HS  = diag([ones(P.Dim,1);zeros(nLamS,1)]);
H   = [HP HPS; HPS' HS];
qp.H = H;
qp.f = [];

sol = mpt_solve(qp);
ret.exitflag = sol.exitflag;

if sol.exitflag==MPTOPTIONS.OK
    y = sol.xopt(1:P.Dim);
    x = sol.xopt(size(matP.A,2)+(1:P.Dim));
    ret.dist = norm(x-y);
    ret.x = x;
    ret.y = y;
end

end
