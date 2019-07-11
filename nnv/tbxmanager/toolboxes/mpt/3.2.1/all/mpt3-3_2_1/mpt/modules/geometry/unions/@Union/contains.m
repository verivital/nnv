function [isin, inwhich, closest] = contains(U, x, fastbreak)
%
% Checks whether the union contains a given point
%
% [isin, inwhich, closest] = U.contains(x, [fastbreak])
%
% inputs:
%   U: single union (use U.forEach() to evaluate arrays), mandatory
%   x: point to check, mandatory
%   fastbreak: fast-abortion boolean flag, false by default
%
% outputs:
%   isin: true if at least one set of the union contains "x"
%   inwhich: indices of sets that contain "x" (will be a singleton if
%            "fastbreak=true")
%   closest: if no set contains "x", this output contains index of the set
%            which is closest to "x". Note: since this computation is
%            expensive, do not ask for the third output unless you really
%            need it.

narginchk(2, 3);
% use obj.forEach(@(u) u.contains(x)) to evaluate arrays
error(U.rejectArray());

if numel(U)==0
	isin = [];
	inwhich = [];
	closest = [];
	return
end
if nargin<3
    fastbreak = false;
end

%% validation
validate_realvector(x);
if size(x, 2)~=1
	error('The point must be a column vector.');
end
nx = numel(x);
iscell_set = iscell(U.Set);
for i = 1:U.Num
	if ( iscell_set && U.Set{i}.Dim ~= nx ) || ...
			( ~iscell_set && U.Set(i).Dim ~= nx )
		error('All sets must be %d-dimensional.', nx);
	end
end

%% search
isin = false;
inwhich = [];
closest = [];
for i = 1:U.Num
	if ( iscell_set && U.Set{i}.contains(x) ) || ...
			(~iscell_set && U.Set(i).contains(x))
		isin = true;
		inwhich = [inwhich, i];
		if fastbreak
			return
		end
	end
end

%% find closest region if necessary
if ~isin && nargout>2
	d = Inf(1, U.Num);
	for i=1:U.Num
		if iscell_set
			s = distance(U.Set{i}, x);
		else
			s = distance(U.Set(i), x);
		end
		if ~isempty(s.dist)
			d(i) = s.dist;
		end
	end
	[~, closest] = min(d);
end

end
