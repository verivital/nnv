function filter = filter_softMax(obj)
% Soft upper bound constraint

global MPTOPTIONS
if isempty(MPTOPTIONS)
	MPTOPTIONS = mptopt;
end

% set up the filter
filter = FilterSetup;
filter.addField('penalty', AffFunction(MPTOPTIONS.infbound*ones(1, obj.n)), @(x) isa(x, 'Function'));
filter.addField('maximalViolation', 1e3*ones(obj.n, 1), @isnumeric);

% this filter depends on the "max" filter
filter.dependsOn('max') = true;

% the filter impacts the following calls:
filter.callback('constraints') = @on_constraints;
filter.callback('objective') = @on_objective;
filter.callback('instantiate') = @on_instantiate;
filter.callback('uninstantiate') = @on_uninstantiate;
filter.callback('addFilter') = @on_addFilter;
filter.callback('removeFilter') = @on_removeFilter;
filter.callback('getVariables') = @on_variables;
filter.callback('set') = @on_set;

end

%------------------------------------------------
function out = on_variables(obj, varargin)
% called when filter's variables are requested

% Response: structure (or an array of structures) with following fields:
%
%  .var: sdpvar representation of the introduced variable
%  .parametric: logical, if true, the variable will become part of the
%              vector of initial conditions
if isa(obj.Internal.soft_max, 'sdpvar')
	out.var = obj.Internal.soft_max;
	out.parametric = false;
else
	out = [];
end

end

%------------------------------------------------
function out = on_instantiate(obj, varargin)
% called after the object was instantiated

% soft constraint require introducing new variables
obj.Internal.soft_max = sdpvar(obj.n, obj.N, 'full');
out = [];

end

%------------------------------------------------
function out = on_uninstantiate(obj, varargin)
% called when the YALMIP representation of variables is removed

% soft constraint require introducing new variables
obj.Internal.soft_max = [];
out = [];

end

%------------------------------------------------
function out = on_constraints(obj, varargin)
% called when creating constraints

s = obj.Internal.soft_max;
out = [];
for i = 1:obj.n
	out = out + [ obj.var(i, :) <= obj.max(i, :) + s(i, :) ];
	out = out + [ 0 <= s(i, :) <= obj.softMax.maximalViolation(i, :) ];
end

end

%------------------------------------------------
function out = on_objective(obj, varargin)
% called when creating the objective function

s = obj.Internal.soft_max;
out = 0;
for k = 1:obj.N
	out = out + obj.softMax.penalty.feval(s(:, k));
end
        
end

%------------------------------------------------
function obj = on_addFilter(obj)
% called after the filter was added

% we need to deactivate the "max" filter
obj.disableFilter('max');

end

%------------------------------------------------
function obj = on_removeFilter(obj)
% called prior to the filter is removed

% we need to re-activate the "max" filter
obj.enableFilter('max');

end

%------------------------------------------------
function obj = on_set(obj, value)
% validation

error(validate_vector(value.maximalViolation, obj.n, 'maximal violation'));
error(obj.validatePenalty(value.penalty));

obj.softMax.maximalViolation = value.maximalViolation;
obj.softMax.penalty = value.penalty;

end
