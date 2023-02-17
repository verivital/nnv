function z = tfAdd(x, y)
import mlp_nn.ops.*;

%   Copyright 2020-2021 The MathWorks, Inc.

xrank = x.rank; 
yrank = y.rank; 
zrank = max(xrank, yrank); 
x = x.value; 
y = y.value; 

if ~isfloat(x) && isfloat(y)
x = cast(x, 'like', y); 
elseif ~isfloat(y) && isfloat(x) 
y = cast(y, 'like', x); 
end 

if isscalar(x) || isscalar(y)
z.value = x + y;
z.rank = zrank; 
return; 
end 
% in MATLAB, broadcasting starts from the left. In TF, broadcasting 
% starts from the right. For this reason, we will permute x and y to 
% match the reverse TF dimension. Using rank, we will know with certainty 
% the correct permutation to reverse dimensions. We will assume non-dlt
% formatted arrays are already in reverse TF dimension. 
isXDLTFormat = isa(x, 'dlarray') && ~isempty(x.dims) && ~all(x.dims == 'U') && xrank > 1; 
isYDLTFormat = isa(y, 'dlarray') && ~isempty(y.dims) && ~all(y.dims == 'U') && yrank > 1;
if isXDLTFormat 
[xPermutationVec, xtflabel] = sortToTFLabel(1:ndims(x), x.dims); 
x = stripdims(x);
x = permute(x, flip(xPermutationVec)); 
elseif isa(x, 'dlarray')
x = stripdims(x); 
end 

if isYDLTFormat
[yPermutationVec, ytflabel] = sortToTFLabel(1:ndims(y), y.dims); 
y = stripdims(y);
y = permute(y, flip(yPermutationVec)); 
elseif isa(y, 'dlarray')
y = stripdims(y); 
end

z = x + y; 
if isXDLTFormat && isYDLTFormat && numel(xtflabel) == numel(ytflabel) && all(xtflabel == ytflabel)
if zrank > 1
z = permute(z, zrank:-1:1); 
end 
z = dlarray(z, xtflabel);
elseif isXDLTFormat && xrank == zrank
if zrank > 1
z = permute(z, zrank:-1:1); 
end 
z = dlarray(z, xtflabel); 
elseif isYDLTFormat && yrank == zrank
if zrank > 1
z = permute(z, zrank:-1:1); 
end 
z = dlarray(z, ytflabel); 
else 
z = dlarray(z, repmat('U', [1 ndims(z)])); 
end
z = struct('value', z, 'rank', zrank); 
end 
