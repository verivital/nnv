function filter = filter_block(varargin)
% Filter's skeleton

% set up the filter
filter = FilterSetup;
filter.addField('from', [], @isnumeric);
filter.addField('to', [], @isnumeric);

% the filter impacts the following calls:
filter.callback('constraints') = @on_constraints;

end

%------------------------------------------------
function out = on_constraints(obj, varargin)
% called when constructing constraints

out = [];
for k = obj.block.from:obj.block.to-1
    out = out + [ obj.var(:, k) == obj.var(:, k+1) ];
end

end
