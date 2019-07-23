function ts = contains(obj,x)
%
% YSet/contains
%
% Synopsis
% ---
%
% Tests whether a point is contained in a YSet object
%
% Syntax
% ---
%
% status = Y.contains(x)
%
% Inputs:
% ---
%
% Y: single YSet or an array of YSets
% x: point or set of points
%
% Note: if "x" is a double, then it must be either a column vector, or a
% matrix composed of column vectors. No automatic transposition of "x" is
% performed!
%
% Outputs:
% ---
%
% if "ny" is the number of elements in "Y" and "nx" the number of points in
% "x", then:
%   status = (ny x nx) matrix of logicals with "status(i, j)=true" iff
%            "Y(i)" contains point "x(:, j)"

global MPTOPTIONS
if isempty(MPTOPTIONS)
    MPTOPTIONS = mptopt;
end

narginchk(2, 2);
validate_realmatrix(x);

if numel(obj)==0
	ts = [];
	return
end

ny = numel(obj);
nx = size(x, 2);
ts = false(ny, nx);
if obj(1).Dim ~= size(x, 1)
	error('The point must have %d rows.', obj(1).Dim);
end

for i = 1:ny
	for j = 1:nx
		% check residuals
		assign(obj(i).vars, x(:, j));
		residual = checkset(obj(i).constraints);
		ts(i, j) = all(residual>-MPTOPTIONS.abs_tol);
	end
end

end
