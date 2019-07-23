function Q = getFacet(P, varargin)
%
% Returns the facet of the polyhedron P for given inequality index. 
%

narginchk(1, 2);

% empty array
if numel(P)==0
    Q = Polyhedron;    
    return
elseif numel(P)>1
	error('This function does not support arrays of polyhedra. Use forEach().');
end

if ~P.irredundantHRep
    error('Polyhedron must be in its minimal representation. Use "minHRep()" to perform the redundancy elimination.');
end

n1 = size(P.H,1);
if nargin==2
	index = varargin{1}(:)';
	validate_indexset(index);
	if numel(index)>n1 || any(index>n1)
		error('Facet index set contains indices out of range.');
	end
else
	% no index provided => return all facets
	index = 1:n1;
end

Q = [];
for i = index
	% take complement constraints and form H-rep
	ic = setdiff(1:n1, i);

	% form new polyhedron
	Q = [Q; Polyhedron('H',P.H(ic,:),'He',[P.H(i,:); P.He])];
end

end
