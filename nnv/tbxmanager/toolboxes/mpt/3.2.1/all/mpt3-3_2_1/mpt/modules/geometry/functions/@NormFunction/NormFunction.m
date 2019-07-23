classdef NormFunction < Function
	%
	% class for representing norm-based functions
	%
	% syntax:
	%   f = NormFunction(1)      : f = norm(x, 1)
	%   f = NormFunction(Inf)    : f = norm(x, Inf)
	%   f = NormFunction(1, Q)   : f = norm(Q*x, 1)
	%   f = NormFunction(Inf, Q) : f = norm(Q*x, Inf)
	%
	% "Q" need not to be square. Function value is always scalar.
	%
	% 2-norms are not supported because they are neither quadratic, nor
	% piecewise linear.
	
	properties(SetAccess=private)
		weight=1; % weight (1 by default)
		type=1; % either 1 or Inf
		D=0; % dimension of the domain
		R=1; % dimension of the range, norms are always scalar-valued
	end
	
	methods
		
		% Constructor
		function obj = NormFunction(type, Q)
			%
			% syntax:
			%   f = NormFunction(1)      : f = norm(x, 1)
			%   f = NormFunction(Inf)    : f = norm(x, Inf)
			%   f = NormFunction(1, Q)   : f = norm(Q*x, 1)
			%   f = NormFunction(Inf, Q) : f = norm(Q*x, Inf)
			%
			% "Q" need not to be square. Function value is always scalar.
			%
			% 2-norms are not supported because they are neither quadratic,
			% nor piecewise linear.
			%
			% Do not use these objects in the user interface. Use
			% OneNormFunction and InfNormFunction objects instead.
			
			if nargin==0
				% when called from derived classes
				return
			end
			
			narginchk(1, 2);
			
			% validation of arguments is done in setters
			obj.type = type;
			if nargin==2
				obj.weight = Q;
			end
			
			obj.Handle = @(x) norm(obj.weight*x, obj.type);
		end
		
		function obj = set.weight(obj, Q)
			% obj.Q setter
			
			if isempty(Q)
				% restore norm(x, type)
				obj.weight = 1;
				obj.D = 0;
			elseif isscalar(Q)
				% restore unrestricted domain
				validate_realmatrix(Q);
				obj.weight = Q;
				obj.D = 0;
			else
				% domain is equal to number of columns of Q
				validate_realmatrix(Q);
				if issparse(Q)
					Q = full(Q);
				end
				obj.weight = Q;
				obj.D = size(Q, 2);
			end
		end
		
		function obj = set.type(obj, type)
			% obj.type setter
			
			if ~isnumeric(type) || ~ismember(type, [1, Inf])
				error('Norm type can only be either 1 or Inf.');
			end
			obj.type = type;
		end
		
		function display(obj)
			
			% TODO: inherit from Function, see how it's done in
			% AbstractController
            if numel(obj)>1
				fprintf('Array of %d norm functions\n',numel(obj));
				return
            end
			
            if numel(obj)==0
                disp('Empty function');
                return
            end

			if obj.D==0
				fprintf('%s-norm function\n', num2str(obj.type));
			else
				fprintf('%s-norm function in R^%i\n', num2str(obj.type), obj.D);
			end
		end
		
		function status = eq(f, g)
			% Returns true if the two functions are identical
			
			% TODO: move common code to Function/eq
			if numel(f)~=numel(g)
				error('Matrix dimensions must agree.');
			elseif numel(f)>1
				% for arrays
				status = false(size(f));
				for i = 1:numel(f)
					status(i) = (f(i)==g(i));
				end
			else
				% for scalars
				
				status = isa(f, 'NormFunction') && ...
					isa(g, 'NormFunction') && ...
					f.D==g.D && ...
					f.type==g.type && ...
					isequal(f.weight, g.weight);
			end
		end
		
		function status = ne(f, g)
			% Returns true if two functions are not identical
			
			% TODO: inherit from Function
			status = ~eq(f, g);
		end
		
	end
	
end
