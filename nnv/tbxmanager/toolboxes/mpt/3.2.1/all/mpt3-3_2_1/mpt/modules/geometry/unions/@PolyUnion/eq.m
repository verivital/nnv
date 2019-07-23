function ts = eq(U1,U2)
%
% test if two unions of polyhedra cover the same set
%  
%  ts = (U1 == U2)
%  ts = U1.eq(U2)
%

if isa(U1, 'Polyhedron')
	U1 = PolyUnion(U1);
end
if isa(U2, 'Polyhedron')
	U2 = PolyUnion(U2);
end
if ~isa(U1,'PolyUnion') || ~isa(U2,'PolyUnion')
    error('All inputs must be PolyUnion objects.');
end

% both polyhedra are empty arrays
if numel(U1)==0 && numel(U2)==0
    ts = true;
    return;
end
% of the polyhedra is empty array
if numel(U1)==0 || numel(U2)==0
    ts = false;
    return;
end

if numel(U2)~=1
    error('Only single "PolyUnion" object can be tested for equality.');
end

% deal with arrays
nU = numel(U1);
if numel(U1)>1
   ts = false(size(U1));
   for i=1:nU
        ts(i) = U1(i).eq(U2);
   end
   return;
end


if U1.Num==0 && U2.Num==0
    ts = true;
    return;
end
if U1.Num==0 || U2.Num==0
    ts = false;
    return
end

% check dimensions
dims = cell(nU, 1);
dims{:} = U1.Dim;

if any(diff([dims{:}]))
    error('The array of polyunions must be in the same dimension.');
end

if U1(1).Dim~=U2.Dim
    error('Unions must have the same dimension.');
end

% use the comparison based on set difference for arrays of polyhedra
ts = (U1.Set == U2.Set);


end
