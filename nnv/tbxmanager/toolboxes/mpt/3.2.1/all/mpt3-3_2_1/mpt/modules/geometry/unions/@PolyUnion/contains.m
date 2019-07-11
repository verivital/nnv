function [isin, inwhich, closest] = contains(obj, x, fastbreak)
%
% Checks whether the union contains a given point
%
% [isin, inwhich, closest] = U.contains(x, [fastbreak])
%
% inputs:
%   U: single polyunion (use U.forEach() to evaluate arrays), mandatory
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
error(obj.rejectArray());

%% parse inputs
if nargin<3
	fastbreak = false;
end
check_hull_first = false;

isin = [];
inwhich = [];
closest = [];
if numel(obj)==0 || ( numel(obj)==1 && obj.Num==0 )
	return
end

%% validation
validate_realvector(x);
if size(x, 2)~=1 || numel(x)~=obj.Dim
	error('The point must be a %dx1 vector.', obj.Dim);
end

%% heuristics
if check_hull_first
	% check the convex hull first, maybe we could exit quickly.
	%
	% in theory, this is a great idea. in practice, it does not pay out.

	if numel(obj)==1 && isfield(obj.Internal, 'convexHull') && ...
			~isempty(obj.Internal.convexHull)
		H = obj.Internal.convexHull;
		if ~H.contains(x)
			% not in the convex hull, no point in evaluating further
			isin = false;
			if nargout==3
				% computing the closest region is slow anyhow, so there is
				% no harm in going via Union/contains
				[~, ~, closest] = obj.contains@Union(x, fastbreak);
			end
			return
		end
	end
end

%% search
% exploit Polyhedron/contains operating on arrays
c = obj.Set.contains(x, fastbreak);
isin = any(c);
% always return "inwhich" as a row vector
inwhich = find(c); 
inwhich = inwhich(:)';

%% find closest region if necessary
if nargout==3 && ~isin
	% computing the closest region is slow anyhow, so there is no harm in
	% going via Union/contains
	[~, ~, closest] = obj.contains@Union(x, true);
end

end
