function filter = filter_penalty(varargin)
% Penalizes the signal in the cost function
%
% If the signal is a state variable, this filter DOES NOT penalize the
% final predicted state (aka the "terminal state").

% set up the filter
filter = FilterSetup;
filter.addField('value', []);

% the filter impacts the following calls:
filter.callback('objective') = @on_objective;
filter.callback('set') = @on_set;

end

%------------------------------------------------
function out = on_objective(obj, varargin)
% called when constructing the cost function

out = 0;
if isempty(obj.penalty)
	return
end

% remember that state variables have length N+1, but we only
% penalize the first N components. the terminal penalty can be
% easily added by obj.with('terminalPenalty')
if obj.isKind('x')
	M = obj.N-1;
else
	M = obj.N;
end

reference = zeros(obj.n, M);
if obj.hasFilter('reference')
	if ismember(obj.Internal.reference.type, {'free', 'preview', 'symbolic'})
		% symbolic reference (vector or matrix)
		reference = obj.Internal.reference.var;
	elseif ~isempty(obj.reference)
		% numerical reference (vector or matrix)
		reference = obj.reference;
	end
end

for k = 1:M
	% if "k" exceeds number of references, repeat with the last provided
	% reference. If size(reference,2)>1, then we get trajectory preview.
	k_ref = min(k, size(reference, 2));
	value = obj.var(:, k) - reference(:, k_ref);
	out = out + obj.penalty.feval(value);
end

end

%------------------------------------------------
function obj = on_set(obj, P)
% called prior to property being set

% validate the penalty (empty penalty means no penalization)
if ~isempty(P)
	error(obj.validatePenalty(P));
end
obj.penalty = P;

end
