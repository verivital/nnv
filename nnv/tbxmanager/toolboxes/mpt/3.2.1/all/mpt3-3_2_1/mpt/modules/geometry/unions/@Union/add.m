function obj = add(obj, C)
%
%  Add the object to the union
%

% deal with arrays
if numel(obj)>1
    for i=1:numel(obj)
        obj(i) = obj(i).add(C);
    end
    return;
end

if ~isa(C,'ConvexSet')
    error('The argument must be derived from "ConvexSet" class.');
end

% remove empty sets
c = isEmptySet(C);
C(c) = [];

% check function hadles
fnames = obj.listFunctions;
nf = numel(fnames);
N = length(C);

% if union contains some function handles, the new added sets must also
% have the same function handles
if numel(fnames)>0 && N>0
	for i = 1:N
		if any(~C(i).hasFunction(fnames))
            error('The set %i to be added holds different function names than the union.',i);
		elseif length(C(i).Functions) ~= nf
			error('All sets to be added must have associated %d number of functions.',nf);
		end
	end
end

if N>0
	obj.Set(obj.Num+1:obj.Num+N,1) = num2cell(C(:));
	%obj.Set(obj.Num+1:obj.Num+N,1) = C(:);
end
   
end
