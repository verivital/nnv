function y = tfMatMul(A, B, transp_a, transp_b)
import mlp_nn.ops.*;

%   Copyright 2020-2021 The MathWorks, Inc.

A = A.value; 
B = B.value; 

if transp_a
A = A';
end

if transp_b
B = B';
end


if isa(A, 'dlarray')
alabels = A.dims;
ndimsa = ndims(A);
A = stripdims(A);
else
alabels = [];
ndimsa = ndims(A);
end

if isa(B, 'dlarray')
blabels = B.dims;
ndimsb = ndims(B); 
B = stripdims(B);
else
blabels = [];
ndimsb = ndims(B); 
end

% non-formatted dlarray data will be assumed to have reverse TF format. 
% So and no permutation will occur. 
if ~isempty(alabels) && all(alabels ~= 'U')
[permutea, tfalabels] = sortToTFLabel(1:ndimsa, alabels); 
permutea = fliplr(permutea);
A = permute(A, permutea); 
end 

if ~isempty(blabels) && all(blabels ~= 'U')
[permuteb, tfblabels] = sortToTFLabel(1:ndimsb, blabels); 
permuteb = fliplr(permuteb);
B = permute(B, permuteb);
end

% Thm: y' = (A * B)' == B' * A'. 
% the input matrices are in reverse TensorFlow format now. (lets call them A',
% and B'). Forward TensorFlow format is called A and B respectively. in TensorFlow, the last
% two dimensions are multiplied. pagemtimes will multiply the first two
% dimensions, Everything else is a page dimension. By the Thm above, we
% will call pagemtimes(B', A'). This will return a transposed y. We then 
% apply the reversed TensorFlow labels to get the correct MATLAB labels. If 
% labels are unknown, we will output the reverse TensorFlow dimension tensor.
y = pagemtimes(B, A);

yrank = ndims(y); 
% re-apply labels to the reverse TensorFlow dimension tensor 
if ~isempty(alabels) && ~isempty(blabels) && ~any(blabels == 'U') && ~any(alabels == 'U')
% Apply back the reverse TF labels. 
tfalabels = fliplr(tfalabels); 
tfblabels = fliplr(tfblabels); 
ylabel = fliplr([tfblabels(1) tfalabels(2:end)]); 
if yrank > 1
y = permute(y, yrank:-1:1); 
end 
y = dlarray(y, ylabel);
elseif ~isempty(alabels) && all(alabels ~= 'U') && ~isempty(tfalabels)
% Apply output labels. 
if yrank > 1
y = permute(y, yrank:-1:1); 
end 
y = dlarray(y, tfalabels); 
elseif ~isempty(blabels) && all(blabels ~= 'U') && ~isempty(tfblabels) 
% Try to re-apply the labels that came from MATLAB. 
if yrank > 1
y = permute(y, yrank:-1:1); 
end 
y = dlarray(y, tfblabels); 
end

y = struct('value', y, 'rank', yrank); 
end 
