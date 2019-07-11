function sep = separate(obj, x)
% SEPARATE Compute a separating hyperplane between this set and the
%          given point x.
%
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
%    a,b  : a is vector of length P.Dim and b is a scalr
%           s.t. a'*x <= b and a'*y >= b for all y in P
%  or [] if x in P or P is empty.
%

global MPTOPTIONS

if isempty(MPTOPTIONS)
    MPTOPTIONS = mptopt;
end

narginchk(2, 2);
error(obj.rejectArray());

validate_realvector(x);

if numel(x)~=obj.Dim
    error('The point "x" must be a vector with the length of %i.', obj.Dim);
end

% for empty objects, quickly return
if obj.isEmptySet
    sep = [];
    return;
end

x = x(:);
ret = obj.project(x);
switch ret.exitflag
    case MPTOPTIONS.INFEASIBLE,
        sep = [];
        return
    case MPTOPTIONS.UNBOUNDED,
        error('Solver returned unbounded. This should not happen');
end
if ret.dist < MPTOPTIONS.rel_tol,
    sep = [];
    return;
end

y = ret.x; % Projected point

sep = y-x; % Hyperplane normal points from x to y
sep(end+1) = sep'*(y+x)/2; % Hyperplane intersects point halfway between x and y
sep = sep(:)';

end
