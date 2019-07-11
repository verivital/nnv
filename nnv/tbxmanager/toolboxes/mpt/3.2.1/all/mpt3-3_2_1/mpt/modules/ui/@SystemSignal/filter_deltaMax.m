function filter = filter_deltaMax(obj)
% Upper bound on the slew rate

% set up the filter
filter = FilterSetup;
filter.addField('value', Inf(obj.n, 1), @isnumeric);

% the filter impacts the following calls:
filter.callback('constraints') = @on_constraints;
filter.callback('addFilter') = @on_addFilter;
filter.callback('removeFilter') = @on_removeFilter;
filter.callback('set') = @on_set;

end

%-----------------------------------------------
function obj = on_set(obj, value)
% called when the filter's values are changed

error(validate_vector(value, obj.n, 'value of deltaMax'));
if obj.hasFilter('deltaMin') && any(value < obj.deltaMin)
	error('Upper bound cannot exceed lower bound.');
end
obj.deltaMax = value;

end

%------------------------------------------------
function out = on_constraints(obj, varargin)
% called when creating constraints

out = [];
if ~obj.isKind('x')
	% non-state signals have the previous value in
	% obj.Internal.previous.var
	previous = obj.Internal.previous.var;

	% constraints
	delta = obj.var(:, 1)-previous;
	for i = 1:obj.n
		% Do not include +/-Inf bounds
		if ~isinf(obj.deltaMax(i))
			out = out + [ delta(i) <= obj.deltaMax(i) ];
		end
	end
end

for k = 2:obj.N
	delta = obj.var(:, k)-obj.var(:, k-1);
	for i = 1:obj.n
		% Do not include +/-Inf bounds
		if ~isinf(obj.deltaMax(i))
			out = out + [ delta(i) <= obj.deltaMax(i) ];
		end
	end
end	

end

%------------------------------------------------
function out = on_addFilter(obj, varargin)
% called after the filter was added

% non-state signals require the previous value to be stored in vector of
% initial conditions
if ~obj.isKind('x')
	if ~obj.hasFilter('previous')
		obj.addFilter('previous');
	end
	% register ourselves with the "previous" filter (only register
	% once)
	obj.previous = union(obj.previous, mfilename);
end
out = [];

end

%------------------------------------------------
function out = on_removeFilter(obj, varargin)
% called prior to the filter is removed

% non-state signals require the previous value to be stored in vector of
% initial conditions
if ~obj.isKind('x')
	% deregister ourselves from the "previous" filter
	obj.previous = setdiff(obj.previous, mfilename);
	if isempty(obj.previous)
		% we were the last ones requiring the previous value, remove the
		% "init" filter
		obj.removeFilter('previous');
	end
end
out = [];

end
