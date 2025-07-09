function [indexingCell, numDimsY] = prepareSliceArgs(X, Starts_, Ends_, Axes__, Steps__, numDimsX_)
% Prepares arguments for implementing the ONNX Slice operator
%#codegen

%   Copyright 2024 The MathWorks, Inc.

% Starts, Ends and Axes are all origin 0. Axes refer to the ONNX dimension
% ordering, but X uses the reverse, DLT ordering. Starts, Ends, Axes, and
% Steps correspond positionally. Axes and Steps may be omitted, with
% defaults described in the ONNX spec.

Starts = relu_quadrotor2d_state.coder.ops.extractIfDlarray(Starts_);
Ends = relu_quadrotor2d_state.coder.ops.extractIfDlarray(Ends_);
Axes_ = relu_quadrotor2d_state.coder.ops.extractIfDlarray(Axes__);
Steps_ = relu_quadrotor2d_state.coder.ops.extractIfDlarray(Steps__);
numDimsX = relu_quadrotor2d_state.coder.ops.extractIfDlarray(numDimsX_);

% Set default Steps and Axes if not supplied
if isempty(Steps_)
    Steps = ones(1, numel(Starts));   % Size change
else
    Steps = Steps_;
end

if isempty(Axes_)
    Axes = 0:numDimsX-1;   % All axes
else
    Axes = Axes_;   
    % Handle negative Axes.
    coder.unroll();
    for i = 1:numel(Axes)
        if Axes(i)<0
            Axes(i) = Axes(i) + numDimsX;
        end
    end
end

% Convert to DLT Axes.
DLTAxes = Axes;
coder.unroll();
for i = 1:numel(Axes)
    DLTDim = numDimsX - Axes(i);                                               % The DLT dim is the reverse of the ONNX dim.
    DLTAxes(i) = DLTDim;  
    % "If a negative value is passed for any of the start or end indices,
    % it represents number of elements before the end of that dimension."
    if Starts(i) < 0
        Starts(i) = size(X,DLTDim) + Starts(i);
    end
    if Ends(i) < 0
        Ends(i) = max(-1, size(X,DLTDim) + Ends(i));                        % The -1 case is when we're slicing backward and want to include 0.
    end
    % "If the value passed to start or end is larger than the n (the number
    % of elements in this dimension), it represents n."
    if Starts(i) > size(X,DLTDim)
        Starts(i) = size(X,DLTDim);
    end
    if Ends(i) > size(X,DLTDim)
        Ends(i) = size(X,DLTDim);
    end
end
indexingCell = cell(1, numDimsX);

coder.unroll();
for i=1:numDimsX
    j = nnet.internal.cnn.coder.onnx.find(DLTAxes == i);
    if ~isempty(j)
        % Use ONNX indexing
        if Steps(j) > 0
            % numelIdx = ceil((Ends(j) - Starts(j))/Steps(j));
            % idx = zeros(1,numelIdx);
		    idx = (Starts(j) : Steps(j) : Ends(j)-1);
            indexingCell{i} = 1 + idx;            % 1 + (Origin 0 indexing with end index excluded)
        else
            % numelIdx = ceil((Ends(j) - Starts(j) +2)/Steps(j));
            % idx = zeros(1,numelIdx);
		    idx = (Starts(j) : Steps(j) : Ends(j)+1);
            indexingCell{i} = 1 + idx;            % 1 + (Origin 0 indexing with end index excluded)
        end
    else
        % Select all elements at dim
        indexingCell{i} = 1:size(X,i);
    end
end
numDimsY = numDimsX;
end
