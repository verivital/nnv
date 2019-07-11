function ts = isFullDim(obj)
%
% check if all polyhedra in the set are full-dimensional
%

% deal with arrays
if numel(obj)>1
    ts = -ones(size(obj));
    for i=1:numel(obj)
        ts(i) = obj(i).isFullDim;
    end
    return;
end

% empty obj
if obj.Num==0
    ts = false;
    return;
end

% check full-dimensionality
if isempty(obj.Internal.FullDim)            
    ts = all(obj.Set.isFullDim);
    obj.Internal.FullDim = ts;
else
   ts = obj.Internal.FullDim; 
end



end
