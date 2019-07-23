classdef InfNormFunction < NormFunction
	%
	% represents weighted inf-norm functions
	%
	% syntax:
	%   f = InfNormFunction(Q) : f = norm(Q*x, Inf)
	%
	% "Q" need not to be square. Function value is always scalar.
	
	methods
		
		% Constructor
		function obj = InfNormFunction(Q)
			% Constructs a weighted 1-norm function object
			%
			% syntax:
			%   f = InfNormFunction(Q) : f = norm(Q*x, Inf)
			%
			% "Q" need not to be square. Function value is always scalar.
			
			narginchk(1, 1);
			obj = obj@NormFunction(Inf, Q);
		end
		
	end
	
end
