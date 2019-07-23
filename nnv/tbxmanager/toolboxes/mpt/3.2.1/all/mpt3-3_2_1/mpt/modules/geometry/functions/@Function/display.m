function display(obj)
%
% overwrites display for Function object
%

if length(obj)>1
    fprintf('Array of %d Functions.\n',numel(obj));
    return
end

if isempty(obj.Handle)
    disp('Empty Function');
else
    s=func2str(obj.Handle);
    l = length(s);    
    if l>50
        c = [s(1:50),'...'];
    else
        c=s;
    end
    disp(['Function: ',c]);
end

end