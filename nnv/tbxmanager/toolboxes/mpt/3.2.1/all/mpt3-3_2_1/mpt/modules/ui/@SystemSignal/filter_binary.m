function filter = filter_binary(obj)
% Binary character of the signal

% set up the filter
filter = FilterSetup;
filter.addField('index', 1:obj.n, @isnumeric);

% the filter impacts the following calls:
filter.callback('constraints') = @on_constraints;
filter.callback('set') = @on_set;

end


%------------------------------------------------
function obj = on_set(obj, new)
% called when the filter's parameters are changes

if islogical(new) && new
	% true = all elements are binary
	obj.binary = 1:obj.n;
elseif islogical(new) && ~new
	% false = no binary elements
	obj.binary = [];
else
	% just set indices
	obj.binary = new;
end

end

%------------------------------------------------
function out = on_constraints(obj, varargin)
% called when constructing constraints

if islogical(obj.binary)
	idx = 1:obj.n;
else
	idx = obj.binary;
end

out = [ binary(obj.var(obj.binary, :)) ];
out = out + [ 0 <= obj.var(obj.binary, :) <= 1 ];

end
