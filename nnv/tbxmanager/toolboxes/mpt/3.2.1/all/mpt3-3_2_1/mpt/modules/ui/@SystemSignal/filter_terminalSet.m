function filter = filter_terminalSet(obj)
% Polyhedral terminal set constraint

% terminal set can only be added on state variables
if ~obj.isKind('x')
	error('Filter "terminalSet" can only be added to state variables.');
end

% set up the filter
filter = FilterSetup;
filter.addField('value', [], @validate_polyhedron);

% the filter impacts the following calls:
filter.callback('constraints') = @on_constraints;
filter.callback('set') = @on_set;

end

%------------------------------------------------
function out = on_constraints(obj, varargin)
% called when constructing constraints

% exit immediately if no terminal set is provided
if ~isa(obj.terminalSet, 'Polyhedron')
	out = [];
	return
end

H = obj.terminalSet.A;
K = obj.terminalSet.b;
Heq = obj.terminalSet.Ae; 
Keq = obj.terminalSet.be;

% TODO: support multiple terminal sets
out = [ H*obj.var(:, end) <= K ];
if ~isempty(Heq)
    out = out + [ Heq*obj.var(:, end) == Keq ];
end

end

%------------------------------------------------
function obj = on_set(obj, P)
% called before the terminal set is changed

% empty penalty means no penalization
if isempty(P)
    return
elseif ~isa(P, 'Polyhedron')
	error('The input must be a polyhedron.');
elseif numel(P)~=1
	error('The input must be a single polyhedron.');
elseif P.Dim ~= obj.n
	error('The polyhedron must be in dimension %d.', obj.n);
elseif P.isEmptySet()
	error('The polyhedron must not be empty.');
end

% we require the H-representation, which we also normalize to avoid
% numerical problems
Q = P.copy(); % make a copy
for i = 1:numel(Q)
	Q(i).minHRep().normalize();
end

obj.terminalSet = Q;

end
