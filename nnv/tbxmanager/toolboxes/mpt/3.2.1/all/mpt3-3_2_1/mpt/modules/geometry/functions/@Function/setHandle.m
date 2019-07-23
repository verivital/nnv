function obj=setHandle(obj, h)
%
% assigns function handle to a "Function" object 
%

narginchk(2, 2);

if numel(h)~=numel(obj)
    error('There must be %d handles present in a cell.',numel(obj));
end

if numel(obj)>1 && ~isa(h,'cell')
    error('Handles must be given in a cell.');
end

if numel(h)==1
    h = {h};
end

% make handles a column vector
h = h(:);

for i=1:numel(obj)
    
    if ~isa(h{i},'function_handle')
        error('Argument must a function handle.')
    end

    
    if strcmp('AffFunction',class(obj(i))) || strcmp('QuadFunction',class(obj(i)))
        error('setHandle: Can not change handle for affine or quadratic functions.');
    end
    
    obj(i).Handle = h{i};
    
end
end
