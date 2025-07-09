function [Y, numDimsY] = onnxConcat(ONNXAxis_, XCell_, numDimsXArray_)
% Concatentation that treats all empties the same. Necessary because
% dlarray.cat does not allow, for example, cat(1, 1x1, 1x0) because the
% second dimension sizes do not match.
%#codegen

% Copyright 2024 The MathWorks, Inc.

    ONNXAxis   = yolov5nano_LRelu_640.coder.ops.extractIfDlarray(ONNXAxis_);
    numDimsXArray   = yolov5nano_LRelu_640.coder.ops.extractIfDlarray(numDimsXArray_);
    
    numDimsY = numDimsXArray(1);
    
    if isempty(XCell_)
        Y = dlarray([]);
    else
        numTensors = coder.const(numel(XCell_));    
        XCell = cell(1, numTensors);
	    coder.unroll();
        for i = 1:numTensors
            if isempty(XCell_{i})
                XCell{i} = [];
            else
                XCell{i} = XCell_{i};
            end
        end
        if ONNXAxis<0
            ONNXAxis = ONNXAxis + numDimsY;
        end
        DLTAxis = coder.const(numDimsY - ONNXAxis);
        Y = cat(DLTAxis, XCell{:});
    end
end
