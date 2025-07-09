function [perm, numDimsA] = prepareTransposeArgs(ONNXPerm_, numDimsA_)
% Prepares arguments for implementing the ONNX Transpose operator
%#codegen

%   Copyright 2024 The MathWorks, Inc.
ONNXPerm = yolov5nano_LRelu_640.coder.ops.extractIfDlarray(ONNXPerm_);
numDimsA = yolov5nano_LRelu_640.coder.ops.extractIfDlarray(numDimsA_);

if numDimsA <= 1        % Tensors of numDims 0 or 1 are unchanged by ONNX Transpose.
    perm = [];
else
    if isempty(ONNXPerm)        % Empty ONNXPerm means reverse the dimensions.
        perm = numDimsA:-1:1;
    else
        perm = numDimsA - flip(ONNXPerm);
    end
end
end
