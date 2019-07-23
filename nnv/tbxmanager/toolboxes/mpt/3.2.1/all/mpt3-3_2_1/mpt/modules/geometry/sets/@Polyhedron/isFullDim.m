function tf = isFullDim(P)
% ISFULLDIM Test if polyhedron is full-dimensional.
%
% ------------------------------------------------------------------
% tf = P.isFullDim
%
% Returns: true if full-dim, false otherwise
%

% Deal with arrays of polyhedra
tf = false(size(P));
for i=1:length(P)
	% use a stored information if possible
	fulldim = P(i).Internal.FullDim;
	if ~isempty(fulldim)
		tf(i) = fulldim;
	else
		% if the polyhedron is empty -> not full dimensional
		if P(i).isEmptySet
			P(i).Internal.FullDim = false;
		else
			% compute interior point only if P(i).Empty property has not been set
			sol = P(i).interiorPoint;
			P(i).Internal.FullDim = sol.isStrict && ~isempty(sol.x);
		end
		tf(i) = P(i).Internal.FullDim;
	end
end
end
