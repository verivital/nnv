function [D, numDimsD] = onnxMatMul(A, B_, numDimsA_, numDimsB_)
% Implements the ONNX MatMul operator.
%#codegen

% Copyright 2024 The MathWorks, Inc.

numDimsA   = relu_quadrotor2d_state.coder.ops.extractIfDlarray(numDimsA_);
numDimsB   = relu_quadrotor2d_state.coder.ops.extractIfDlarray(numDimsB_);

% If B is 1-D, temporarily extend it to a row vector
if numDimsB==1
    B = B_(:)';
else
    B = B_;
end
maxNumDims = coder.const(max(numDimsA, numDimsB));
numDimsD = maxNumDims;

if maxNumDims > 2    
    D = pagemtimes(B, A);    
else
    D_ = B * A;
    if numDimsA==1 || numDimsB==1
        D = D_(:);
        numDimsD = 1;
    else
        D = D_;
    end
end
end
