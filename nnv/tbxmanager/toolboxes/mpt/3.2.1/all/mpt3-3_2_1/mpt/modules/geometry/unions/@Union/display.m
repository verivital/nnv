function display(obj)
%
% display function
%

if numel(obj) > 1
    fprintf('Array of %i Unions.\n', numel(obj));
    return
elseif numel(obj)==0
    fprintf('Empty array of Unions.\n');
    return;
end    
if numel(obj.Set)==0
    fprintf('Empty Union.\n');
    return;
end

fprintf('Union of %d convex sets.\n', obj.Num);

% display functions (implemented in Union/displayFunctions)
obj.displayFunctions;

end
