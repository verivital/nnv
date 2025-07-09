function [Y, numDimsY] = onnxGather(X_, ONNXIdx_, ONNXAxis_, numDimsX_, numDimsIdx_)
% Function implementing the ONNX Gather operator
%#codegen

% In ONNX, 'Gather' first indexes into dimension ONNXAxis of data, using
% the contents of ONNXIdx as the indices. Then, it reshapes the ONNXAxis
% into the shape of ONNXIdx.
%   Example 1:
% Suppose data has shape [2 3 4 5], ONNXIdx has shape [6 7], and axis=1.
% The result has shape [2 6 7 4 5].
%   Example 2:
% Suppose data has shape [2 3 4 5], ONNXIdx has shape [6], and axis=1.
% The result has shape [2 6 4 5].
%   Example 3:
% Suppose data has shape [2 3 4 5], ONNXIdx has shape [] (a scalar), and axis=1.
% The result has shape [2 4 5].
%
% Since we're using reverse indexing relative to ONNX, in this function
% data and ONNXIdx both have reversed dimension ordering.

% Copyright 2024 The MathWorks, Inc.
    
    X          = relu_quadrotor2d_state.coder.ops.extractIfDlarray(X_);
    ONNXIdx    = relu_quadrotor2d_state.coder.ops.extractIfDlarray(ONNXIdx_);
    ONNXAxis   = relu_quadrotor2d_state.coder.ops.extractIfDlarray(ONNXAxis_);
    numDimsX   = relu_quadrotor2d_state.coder.ops.extractIfDlarray(numDimsX_);
    numDimsIdx = relu_quadrotor2d_state.coder.ops.extractIfDlarray(numDimsIdx_);

    numDimsY = numDimsIdx + (numDimsX - 1);
    if isempty(X)
        Y = X;
        return;
    end
    % (1) First, do the indexing part of Gather
    if ONNXAxis<0
        ONNXAxis = ONNXAxis + numDimsX;                                 % Axis can be negative. Convert it to its positive equivalent.
    end
    dltAxis = numDimsX - ONNXAxis;                                      % Convert axis to DLT. ONNXAxis is origin 0 and we index from the end
    % ONNXIdx can have negative components. Make them positive.
    sz = size(X, dltAxis);

    coder.unroll();
    for i = 1:numel(ONNXIdx)
        if ONNXIdx(i)<0
            ONNXIdx(i) = ONNXIdx(i) + sz;
        end
    end
    dltIdx  = ONNXIdx + 1;                                 % ONNXIdx is origin-0 in ONNX, so add 1 to get dltIdx

    % Index into data
    sizeX = size(X, 1:numDimsX);
    indexingCell = iGetIndexingCell(numDimsX,sizeX,dltAxis,dltIdx);
    Y_ = X(indexingCell{:});

    % (2) Now do the reshaping part of Gather.
    shapeY = coder.const(size(Y_));
    shape_ = zeros(1, numDimsX);

    if numel(shapeY) > numDimsX        
        coder.unroll();
        for k = 1:numDimsX
            shape_(k) = shapeY(k);
        end
    elseif numel(shapeY) < numDimsX
        coder.unroll();
        for k = 1:numDimsX
            if k <= numel(shapeY)
                shape_(k) = shapeY(k);
            else
                % Compensate for dropped singletons
                shape_(k) = 1;
            end
        end
    else
        shape_ = shapeY;
    end

    if numDimsIdx == 0
        % Delete the indexed dimension
        shape__ = zeros(1, numel(shape_)-1);
        m = 1;
        coder.unroll();
        for l = 1:numel(shape_)
            if l ~= dltAxis
                shape__(m) = shape_(l);
                m = m+1;            
            end
        end
    elseif numDimsIdx > 1
        % Reshape the indexed dimension into the shape of ONNXIdx
        shape__ = [shape_(1:dltAxis-1) size(ONNXIdx, 1:numDimsIdx) shape_(dltAxis+1:end)];
    else
        shape__ = shape_;
    end
    % Extend the shape to 2D so it's valid MATLAB
    if numel(shape__) < 2
        shape = [shape__ ones(1,2-numel(shape__))];
    else
        shape = shape__;
    end
    Y = reshape(Y_, shape);

    function indexingCell = iGetIndexingCell(numDimsX,sizeX,dltAxis,dltIdx)
        indexingCell = cell(1, numDimsX);
        
        coder.unroll();
        for j=1:numDimsX
            if j==dltAxis
                indexingCell{j} = dltIdx(:);
            else
                indexingCell{j} = 1:sizeX(j);             % Same effect as ":"
            end
        end
    end
end