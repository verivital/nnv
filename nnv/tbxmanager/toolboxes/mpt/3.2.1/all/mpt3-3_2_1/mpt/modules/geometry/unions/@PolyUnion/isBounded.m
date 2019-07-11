function ts = isBounded(obj)
%
% check if all polyhedra in the set are bounded
%

% deal with arrays
if numel(obj)>1
    ts = -ones(size(obj));
    for i=1:numel(obj)
        ts(i) = obj(i).isBounded;
    end
    return;
end

% empty obj
if obj.Num==0
    ts = false;
    return;
end

% check boundedness
if isempty(obj.Internal.Bounded)            
    ts = all(obj.Set.isBounded);
    obj.Internal.Bounded = ts;
else
   ts = obj.Internal.Bounded; 
end



end
