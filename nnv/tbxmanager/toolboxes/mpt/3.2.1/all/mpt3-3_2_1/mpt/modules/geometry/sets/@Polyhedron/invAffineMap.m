function Pnew = invAffineMap(P, T, t)
% Compute the inverse image of the polyhedron under the affine map T*x+t
%
%   M = { x | T*x+t \in P }

validate_realmatrix(T);
if nargin<3
	t = zeros(P(1).Dim, 1);
else
    validate_realvector(t);
end

% deal with arrays
if numel(P)>1
	Pnew(size(P)) = Polyhedron;
    for i=1:numel(P)
        Pnew(i) = P(i).invAffineMap(T, t);
    end
    return
end


% TODO: support non-square mappings
if size(T, 1)~=size(T, 2)
	error('Only square mappings supported.');
end

% TODO: deal with the V-representation directly
if ~P.hasHRep
	% we require the H-representation
	P.minHRep();
end

if isempty(P.He_int)
	% faster call if we have no equalities
	Pnew = Polyhedron(P.A*T, P.b-P.A*t);
else
	Pnew = Polyhedron('A', P.A*T, 'b', P.b-P.A*t, ...
		'Ae', P.Ae*T, 'be', P.be - P.Ae*t);
end

end
