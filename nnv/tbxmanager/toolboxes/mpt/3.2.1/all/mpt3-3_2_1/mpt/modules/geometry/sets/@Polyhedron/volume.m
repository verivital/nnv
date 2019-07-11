function vol = volume(P)
%
% Compute the volume of this polyhedron.<p>
% 
% Volume is infinite if the polyhedron is unbounded and zero if it's not full-dimensional.
% 
% @return volume
%
  
% deal with arrays
no = numel(P);
if no>1
    vol = zeros(size(P));
    for i=1:no
        vol(i) = P(i).volume;
    end
    return;
end

% check emptyness and full-dimensionality
if P.isEmptySet || ~P.isFullDim
    vol = 0;
    return;
end

% computing volume requires vertices
P.minVRep();

if ~P.isBounded
	% unbounded polyhedra have infinite volume
	vol = Inf;
	
elseif P.Dim==1
	% issue #71: we need to handle 1D cases manually
	vol = max(P.V) - min(P.V);

elseif size(P.V, 1) == P.Dim+1
	% cheaper volume computation if the set is a fully-dimensional simplex
	% https://en.wikipedia.org/wiki/Simplex#Volume
	S = P.V';
	D = zeros(size(S, 1), size(S, 2)-1);
	for j = 2:size(S, 2)
		D(:, j-1) = S(:, j)-S(:, 1);
	end
	vol = 1/factorial(size(S, 2)-1)*abs(det(D));

else
	% general computation
	[~, vol] = convhulln(P.V);

end
