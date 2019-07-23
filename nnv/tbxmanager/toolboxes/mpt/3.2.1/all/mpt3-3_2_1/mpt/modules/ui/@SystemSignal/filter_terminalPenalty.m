function filter = filter_terminalPenalty(obj)
% Penalty on the terminal state

% terminal penalty can only be added on state variables
if ~obj.isKind('x')
	error('Filter "terminalPenalty" can only be added to state variables.');
end

% set up the filter
filter = FilterSetup;
filter.addField('value', []);

% the filter impacts the following calls:
filter.callback('objective') = @on_objective;
filter.callback('set') = @on_set;

end

%------------------------------------------------
function out = on_objective(obj, varargin)
% called when constructing the objective function

out = 0;
if isempty(obj.terminalPenalty)
	return
end

% no reference by default
reference = zeros(obj.n, 1);

if obj.hasFilter('reference')
	% reference can either be free (sdpvar) or fixed (last column)
	if ismember(obj.Internal.reference.type, {'free', 'symbolic'})
		% symbolic reference, we implicitly assume it's a vector
		reference = obj.Internal.reference.var;
	elseif ~isempty(obj.reference)
		% fixed reference
		reference = obj.reference(:, end);
	end
end

out = obj.terminalPenalty.feval(obj.var(:, end) - reference);

end

%------------------------------------------------
function obj = on_set(obj, P)
% called prior to property being set

% validate the penalty (empty penalty means no penalization)
if ~isempty(P)
	error(obj.validatePenalty(P));
end
obj.terminalPenalty = P;

end
