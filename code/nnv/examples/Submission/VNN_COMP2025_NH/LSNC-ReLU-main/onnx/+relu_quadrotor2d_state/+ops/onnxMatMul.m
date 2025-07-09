function [D, numDimsD] = onnxMatMul(A, B, numDimsA, numDimsB)
% Implements the ONNX MatMul operator.

% Copyright 2020-2023 The MathWorks, Inc.

% If B is 1-D, temporarily extend it to a row vector
if numDimsB==1
    B = B(:)';
end
maxNumDims = max(numDimsA, numDimsB);
numDimsD = maxNumDims;
if maxNumDims > 2
    % Removes dlarray formats if only one of the input dlarrays is formatted.
    if sum([isempty(dims(A)), isempty(dims(B))]) == 1
        D = pagemtimes(stripdims(B), stripdims(A));
    else
    %computes matrix product of corresponding pages of input arrays A and
    %B.
    D = pagemtimes(B, A);
    end
else
    D = B * A;
    if numDimsA==1 || numDimsB==1
        D = D(:);
        numDimsD = 1;
    end
end
end
