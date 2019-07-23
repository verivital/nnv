function filter = filter_reference(obj)
% Allows the signal to track a certain reference profile

% set up the filter
filter = FilterSetup;
filter.addField('value', zeros(obj.n, 1));

% the filter impacts the following calls:
filter.callback('set') = @on_set;
filter.callback('getVariables') = @on_variables;
filter.callback('instantiate') = @on_instantiate;
filter.callback('uninstantiate') = @on_uninstantiate;
filter.callback('constraints') = @on_constraints;

end

%------------------------------------------------
function out = on_variables(obj, varargin)
% called when filter's variables are requested
%
% Response: structure (or an array of structures) with following fields:
%
%  .var: sdpvar representation of the introduced variable
%  .parametric: logical, if true, the variable will become part of the
%              vector of initial conditions
switch obj.Internal.reference.type
	case {'free', 'preview'}
		out.var = obj.Internal.reference.var(:);
		out.parametric = true;
	case 'symbolic'
		% symbolic references are only used internally (e.g. in
		% LTISystem/filter_integrator), but are not exposed in the vector
		% of initial conditions
		out.var = obj.Internal.reference.var;
		out.parametric = false;
	otherwise
		out = [];
end

end

%------------------------------------------------
function out = on_constraints(obj, varargin)
% called when constructing constraints

out = [];
switch obj.Internal.reference.type
    case 'free'
        % bound symbolic references using signal's min/max values
        %
        % do not include +/-Inf bounds
        v = obj.Internal.reference.var;
        for i = 1:length(v)
            if ~isinf(obj.min(i))
                out = out + [ obj.min(i) <= v(i) ];
            end
            if ~isinf(obj.max(i))
                out = out + [ v(i) <= obj.max(i) ];
            end
        end
    case 'preview'
        % bound symbolic references using signal's min/max values
        %
        % do not include +/-Inf bounds
        v = obj.Internal.reference.var;
        for k = 1:size(v, 2)
            for i = 1:size(v, 1)
                if ~isinf(obj.min(i))
                    out = out + [ obj.min(i) <= v(i, k) ];
                end
                if ~isinf(obj.max(i))
                    out = out + [ v(i, k) <= obj.max(i) ];
                end
            end
        end
end

end

%------------------------------------------------
function out = on_instantiate(obj, varargin)
% called after the object was instantiated

switch obj.Internal.reference.type
    case {'free', 'symbolic'}
        obj.Internal.reference.var = sdpvar(obj.n, 1);
    case 'preview'
        obj.Internal.reference.var = sdpvar(obj.n, obj.N, 'full');
end
out = [];

end

%------------------------------------------------
function out = on_uninstantiate(obj, varargin)
% called when the YALMIP representation of variables is removed

obj.Internal.reference.var = [];
out = [];

end

%------------------------------------------------
function obj = on_set(obj, value)
% called when the reference is to be changed


if isa(value, 'double')
	if ~isempty(value) && (size(value, 1) ~= obj.n)
		error('The refence must have %d rows.', obj.n);
	end
	obj.Internal.reference.type = 'fixed';
	
elseif isa(value, 'char')
	value = lower(value);
	switch lower(value)
		case 'free',
			% free reference which becomes part of the vector of initial
			% conditions
			obj.Internal.reference.type = 'free';
        case 'preview'
            % free reference with trajectory preview becomes part of the
            % vector of initial conditions
            obj.Internal.reference.type = 'preview';
		case 'symbolic',
			% symbolic reference, which is not part of the initial
			% conditions
			obj.Internal.reference.type = 'symbolic';
		otherwise
			error('Unrecognized settings. Can only use "free" or "symbolic".');
	end
	
else
	error('Value of "reference" must be either a double or a string.');
end

obj.reference = value;

end
