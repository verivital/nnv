function filter = filter_min(obj)
% Lower bound constraint

% set up the filter
filter = FilterSetup;
filter.addField('value', -Inf(obj.n, 1), @isnumeric);

% the filter impacts the following calls:
filter.callback('constraints') = @on_constraints;
filter.callback('set') = @on_set;

end

%------------------------------------------------
function out = on_constraints(obj, varargin)
% called when constructing constraints

out = [];
for i = 1:obj.n
	% Do not include +/-Inf bounds
	if any(~isinf(obj.min(i, :)))
		out = out + [ obj.min(i, :) <= obj.var(i, :) ];
	end
end

end

%------------------------------------------------
function obj = on_set(obj, value)
% called when the filter's values are changed

value = value(:);
if numel(value)~=obj.n
	error('Value must be a %dx1 vector.', obj.n);
end
if obj.hasFilter('max') && any(value > obj.max)
	error('Lower bound cannot exceed upper bound.');
end
obj.min = value;

end
