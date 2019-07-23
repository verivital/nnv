function display(obj)
%
% overwrites display for QuadFunction object
%

if numel(obj)>1
    fprintf('Array of %d Quadratic Functions\n',numel(obj));
    return;
end

fprintf('Quadratic Function: R^%d -> R^%d\n',obj.D,obj.R);
%fprintf('y = 0.5*x''*H*x + F*x + g\n');

end