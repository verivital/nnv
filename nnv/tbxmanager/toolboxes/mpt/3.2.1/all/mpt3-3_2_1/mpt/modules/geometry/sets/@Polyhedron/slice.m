function S = slice(P, dims, values, varargin)
%
% Slice a polyhedron by through given dimensions
%
% S = P.slice(dims)
% S = P.slice(dims, values)
% S = P.slice(dims, values, 'keepDim', true/false)
%
% If "values" is omitted, then "values = zeros(size(dims))".
% If "keepDim" is not specified, then "keepDim = false".
%
% Input:
%   P: H-rep polyhedron P = {x | A*x <= b, Aeq*x == beq}
%
% Output:
%   if "keepDim=false":
%     S = { x | A(:, keep)*x <= b - A(:, dims)*values }
%     note that dimension of S is equal to P.Dim-length(dims)
%
%   if "keepDim=true":
%     S = { x | A*x <= b, Aeq*x == beq, x(dims) == values }
%     note that dimension of S is equal to dimension of P.

global MPTOPTIONS

narginchk(2, Inf);
if nargin<3
	values = zeros(size(dims));
end
if nargin<4
	options.keepDim = false;
end
if nargin>3
	% parse options
	ip = inputParser;
	ip.addParamValue('keepDim', false, @islogical);
	ip.parse(varargin{:});
	options = ip.Results;
end


%% deal with arrays
if numel(P)>1
	S = P.forEach(@(e) e.slice(dims, values, varargin{:}));
	return
end
        
%% validation
if P.isEmptySet
    error('Cannot slice empty polyhedra.');
end
% check dimensions
for i=1:numel(dims)
    validate_dimension(dims(i));
end
if any(dims>P.Dim)
    error('The second input cannot exceed dimension of the polyhedron.');
end
if numel(values)~=numel(dims)
	error('"values" must be a vector with %d element(s).', numel(dims));
end

%% computation

% require the H-representation (the getters computes it automatically if
% it does not exist)
if options.keepDim
	% keep dimension of the slice equal to the dimension of the polyhedron
	Z = zeros(numel(dims), P.Dim);
	for i = 1:numel(dims)
		Z(i, dims(i)) = 1;
	end
	S = Polyhedron('H', P.H, 'He', [P.He; [Z, values(:)]]);

else
	% return polyhedron in lower dimension
	keep_dims = setdiff(1:P.Dim, dims);
	remove_dims = dims;
	if isempty(P.He)
		% faster constructor
		S = Polyhedron(P.A(:, keep_dims), P.b-P.A(:, remove_dims)*values(:));
		
	else
		% check that the slice values satisfy equality constraints:
		%
		% find z s.t. P.Ae(:, keep)*z + P.Ae(:, remove)*values == P.be
		nz = length(keep_dims);
		lp.f = zeros(1, nz);
		lp.A = P.A(:, keep_dims);
		lp.b = P.b - P.A(:, remove_dims)*values(:);
		lp.Ae = P.Ae(:, keep_dims);
		lp.be = P.be - P.Ae(:, remove_dims)*values(:);
		lp.lb = []; lp.ub = [];
		lp.quicklp = true;
		sol = mpt_solve(lp);
		
		if sol.exitflag==MPTOPTIONS.OK
			% equality constraints are consistent
			S = Polyhedron('A', P.A(:, keep_dims), ...
				'b', P.b-P.A(:, remove_dims)*values(:), ...
				'Ae', P.Ae(:, keep_dims), ...
				'be', P.be-P.Ae(:, remove_dims)*values(:));
		else
			% equality constraints inconsistent => empty polyhedron of
			% correct dimension
			S = Polyhedron.emptySet(length(keep_dims));
		end
	end
	
	% slice functions
	for n = P.Functions.keys
		name = n{1};
		new = P.Functions(name).slice(dims, values);
		S.addFunction(new, name);
	end
end

end
