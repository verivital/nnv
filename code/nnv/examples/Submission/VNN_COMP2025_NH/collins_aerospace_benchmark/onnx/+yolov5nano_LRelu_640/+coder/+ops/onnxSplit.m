function varargout = onnxSplit(X, ONNXaxis_, splits__, numSplits_, numDimsX_)
% Implements the ONNX Split operator
%#codegen

% Copyright 2024 The MathWorks, Inc.

    ONNXaxis   = yolov5nano_LRelu_640.coder.ops.extractIfDlarray(ONNXaxis_);
    splits_     = yolov5nano_LRelu_640.coder.ops.extractIfDlarray(splits__);
    numSplits  = yolov5nano_LRelu_640.coder.ops.extractIfDlarray(numSplits_);
    numDimsX   = yolov5nano_LRelu_640.coder.ops.extractIfDlarray(numDimsX_);
    
    % ONNXaxis is origin 0. splits is a vector of the lengths of each segment.
    % If numSplits is nonzero, instead split into segments of equal length.
    if ONNXaxis<0
        ONNXaxis = ONNXaxis + numDimsX;
    end
    DLTAxis = coder.const(numDimsX - ONNXaxis);
    
    if numSplits > 0
        C       = size(X, DLTAxis);
        sz      = floor(C/numSplits);
        splits	= repmat(sz, 1, numSplits);
    else
        splits = splits_;
    end
    
    splitIndices = [0 cumsum(splits(:)')];
    sizeX = size(X, 1:numDimsX);
    numY = numel(splitIndices) - 1;
    
    coder.unroll();
    for i = 1:numY
        indexingCell    = cell(1, numDimsX); % Allocate here to avoid overwriting cells with new sizes.
        from            = splitIndices(i) + 1;
        to              = splitIndices(i+1);
        coder.unroll();
        for j=1:numDimsX
            if j==DLTAxis
                indexingCell{j} = from:to;
            else
                indexingCell{j} = 1:sizeX(j);             % Same effect as ":"
            end
        end    
        % The first numY outputs are the Y's. The second numY outputs are their
        % numDims. We assume all the outputs of Split have the same numDims as
        % the input.
        varargout{i}        = X(indexingCell{:});
        varargout{i + numY} = numDimsX;  
    end
end