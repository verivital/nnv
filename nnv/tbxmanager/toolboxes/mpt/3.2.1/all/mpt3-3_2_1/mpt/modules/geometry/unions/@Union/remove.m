function obj = remove(obj, index)
%
%  Remove the object from the union based on the index
%

validate_indexset(index);

% deal with arrays
if numel(obj)>1
    for i=1:numel(obj)
        obj(i) = obj(i).remove(index);
    end
    return;
end

if any(index>obj.Num)
    error('Index is greater than the number of elements in the Union.');
end

% remove based on the index
obj.Set(index) = [];
   
end
