function [DLTShape, numDimsY] = prepareReshapeArgs(X, ONNXShape_, numDimsX, allowzero_)
% Prepares arguments for implementing the ONNX Reshape operator
%#codegen

%   Copyright 2024 The MathWorks, Inc.  

    ONNXShape__ = cGAN_imgSz32_nCh_3.coder.ops.extractIfDlarray(ONNXShape_);
    allowzero   = coder.const(cGAN_imgSz32_nCh_3.coder.ops.extractIfDlarray(allowzero_));
    
    N = coder.const(max(2, numel(ONNXShape__)));
    ONNXShape = flip(ONNXShape__);            % First flip the shape to make it correspond to the dimensions of X.
    
   % Replace zeros with the actual size if allowzero is false
    if allowzero==0 && any(ONNXShape==0)
        i0 = find(ONNXShape==0);
        if numDimsX >= numel(ONNXShape)
            ONNXShape(i0) = size(X, numDimsX - numel(ONNXShape) + i0); % Need to right-align the shape vector and dims
        else
            ONNXShape(i0) = size(X, i0); % Dont need to right-align the shape vector and dims
        end
    end
    
    if isscalar(ONNXShape)
        ONNXShape_ = [ONNXShape 1];
    else
        ONNXShape_ = ONNXShape;
    end
    
    DLTShape = getFixedNewShape(size(X), ONNXShape_, N);
    numDimsY = numel(ONNXShape__);
end

function newShape = getFixedNewShape(curShape, newShape2, outrank)
    newShape = cell(1,numel(newShape2));
    coder.unroll();
    for i = 1:outrank
        if newShape2(i) == -1  
            numKnownSize = prod(newShape2(newShape2~=-1),"all");
            numCurSize = prod(curShape,"all");
            unkSize = numCurSize / numKnownSize;
            newShape{i} = uint32(unkSize);
        else
            newShape{i} = uint32(newShape2(i));
        end 
    end 
end