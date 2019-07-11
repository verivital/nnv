function y=horzcat(varargin)
% Horizontal concatenation of ConvexSet objects
%
% Note: all ConvexSet arrays will be converted to column arrays, regardless
% of their original dimension.

y = vertcat(varargin{:});

end
