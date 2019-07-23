function ts = isConvex(obj)
%
% check if the polyhedron array forms a convex union
%

% deal with arrays
if numel(obj)>1
    ts = -ones(size(obj));
    for i=1:numel(obj)
        ts(i) = obj(i).isConvex;
    end
    return;
end

% empty obj
if obj.Num==0
    ts = false;
    return;
end

if isempty(obj.Internal.Convex)
    % compute the convex hull
	H = obj.convexHull();
	% the union is convex if H\U = emptyset
	ts = all(isEmptySet(mldivide(H, obj.Set, true)));
    
    % store internally
    obj.Internal.Convex = ts;

else
    ts = obj.Internal.Convex;
end



end
