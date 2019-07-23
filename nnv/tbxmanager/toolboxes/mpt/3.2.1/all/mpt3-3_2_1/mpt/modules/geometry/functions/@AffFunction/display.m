function display(obj)
%
% overwrites display for AffFunction object
%
if numel(obj)>1
    fprintf('Array of %d Affine Functions\n',numel(obj));
    return;
end

fprintf('Affine Function: R^%d -> R^%d\n',obj.D,obj.R);
%fprintf('y = F*x + g\n');

end