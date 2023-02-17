function varargout = tfIdentityN(varargin)
import mlp_nn.ops.*;

%   Copyright 2020-2021 The MathWorks, Inc.

varargout = cell(nargin, 1); 
for i = 1:nargin
varargout{i} = varargin{i}; 
end 
end 
