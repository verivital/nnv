function ts = isEmptyFunction(obj)
%
%  if the object contains no function handle, then it is considered as
%  empty
%

ts = zeros(size(obj));
for i=1:numel(obj)
    ts(i) = isempty(obj(i).Handle);
end
