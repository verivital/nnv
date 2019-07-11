function sep = separate(P, S)
% SEPARATE Compute a separating hyperplane between this set and the
%          given set S or point x.
%
% -------------------------------------------------------------------
% sep = P.separate(x)
%
% Projects x onto P, and then builds a hyperplane that separates the
% point x and its projection.
%
% Parameters:
%  x - Vector of length P.Dim
%
% Return:
%  sep - Structure containing the following elements:
%    a,b  : a is vector of length P.Dim and b is a scalar
%           s.t. a'*x <= b and a'*y >= b for all y in P
%  or [] if x in P or P is empty.
%
% -------------------------------------------------------------------
% sep = P.separate(S)
%
% Computes the closest point x in S to P, and then separates x from P.
%
% Parameters:
%  S - Polyhedron
%
% Return:
%  sep - Structure containing the following elements:
%    a,b  : a is vector of length P.Dim and b is a scalar
%           s.t. a'*x <= b for all x in S and a'*y >= b for all y in P
%  or [] if x in P or P is empty.
%

global MPTOPTIONS
if isempty(MPTOPTIONS)
    MPTOPTIONS = mptopt;
end

if isnumeric(S)
    sep = P.separate@ConvexSet(S);
    return;
end

if ~isa(S, 'Polyhedron'),
    error('S must be a polyhedron.');
end
if numel(S)>1
    error('Only single polyhedron S is allowed.');
end

% deal with arrays
error(P.rejectArray());

if S.Dim ~= P.Dim,
    error('S must be in the dimension %i, the same as P.', P.Dim);
end

% Compute the closest points between P and S
ret = distance(P, S);

if ret.exitflag~=MPTOPTIONS.OK
    % infeasible
    sep = [];
    return;
end

% if the sets intersect
if ret.dist < MPTOPTIONS.rel_tol,
    sep = [];
    return;
end


% Compute a separating hyperplane
sep = ret.y-ret.x; % Hyperplane normal points from x to y
sep(end+1) = sep'*(ret.y+ret.x)/2; % Hyperplane intersects point halfway between x and y
sep = sep(:)';


end
