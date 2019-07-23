function supp = support(obj, x)
%
%
% SUPPORT Compute the support of this set in the direction x.
%
% supp = support(x)
%
% Solves the optimization problem:
%
%   support(x) := max x'*y s.t. y in Set
%
% Paramaters:
%  x  - Vector of length "dim" or a matrix with "dim" rows
%
% Returns
%  supp - support of x or [] if empty
%


narginchk(2, 2);

% x can be an array of points put in a matrix
validate_realmatrix(x);

% deal with arrays
if numel(obj)>1
	supp = obj.forEach(@(e) e.support(x));
	return
end

if ~isequal(size(x, 1), obj.Dim)
	error('Input argument "x" must have %d rows.', obj.Dim);
end

% compute the support for each point (points are assumed to be stored
% column-wise)
n_points = size(x, 2);
supp = Inf*ones(n_points, 1);
for i=1:n_points
    sol = obj.extreme(x(:,i));
    if ~isempty(sol.supp)
        supp(i) = sol.supp;
    end
end

end
