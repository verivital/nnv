function dims = prepareReduceArgs(ONNXAxes__, numDimsX_)
% Prepares arguments for implementing the ONNX Reduce operator
%   Copyright 2024 The MathWorks, Inc.
%#codegen

ONNXAxes_ = relu_quadrotor2d_state.coder.ops.extractIfDlarray(ONNXAxes__);
numDimsX = relu_quadrotor2d_state.coder.ops.extractIfDlarray(numDimsX_);

if isempty(ONNXAxes_)
    ONNXAxes = 0:numDimsX-1;   % All axes
else
    ONNXAxes = ONNXAxes_;
    coder.unroll();
    for i = 1:numel(ONNXAxes)
        if ONNXAxes(i)<0
            ONNXAxes(i) = ONNXAxes(i) + numDimsX;
        end
    end
end
dims = numDimsX - ONNXAxes;
end
