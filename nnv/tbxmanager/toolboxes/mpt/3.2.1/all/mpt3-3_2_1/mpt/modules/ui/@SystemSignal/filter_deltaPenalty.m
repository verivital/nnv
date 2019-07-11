function filter = filter_deltaPenalty(varargin)
% Filter's skeleton

% set up the filter
filter = FilterSetup;
filter.addField('value', []);

% the filter impacts the following calls:
filter.callback('objective') = @on_objective;
filter.callback('set') = @on_set;
filter.callback('addFilter') = @on_addFilter;
filter.callback('removeFilter') = @on_removeFilter;

end

%------------------------------------------------
function out = on_addFilter(obj, varargin)

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

%------------------------------------------------
function out = on_objective(obj, varargin)

out = 0;
if isempty(obj.deltaPenalty)
	return
end

if ~obj.isKind('x')
	% non-state signals require a symbolic previous value
	previous = obj.Internal.previous.var;
	% penalization of the first step
	out = out + obj.deltaPenalty.feval(obj.var(:, 1)-previous);
	% now penalize remaining steps
end

% penalize increments
for k = 2:obj.N
	out = out + obj.deltaPenalty.feval(obj.var(:, k) - obj.var(:, k-1));
end

end

%------------------------------------------------
function obj = on_set(obj, P)
% called prior to property being set

% validate the penalty (empty penalty means no penalization)
if ~isempty(P)
	error(obj.validatePenalty(P));
end
obj.deltaPenalty = P;

end
