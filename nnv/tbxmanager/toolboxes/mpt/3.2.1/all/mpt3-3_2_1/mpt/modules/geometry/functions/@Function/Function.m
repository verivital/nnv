classdef Function < handle & IterableBehavior
	%
	% class for representing functions
	%
	% syntax: F = Function('Handle',@fun,'Data',any_data)
	%         F = Function(@fun)
	%
	%         F = Function('Data',struct('p',2));
	%         F.setHandle(@(x)*F.Data.p)
	
	
	%%
	properties(SetAccess = protected)
		Handle % function handle
		Internal % internal data
	end
	
	properties(SetAccess = public)
		Data; % additional user data
	end
	
	%%
	methods(Access = public)
		
		% Constructor
		function F = Function(Handle, Data)
			%
			% sets data for Function object
			%
			% syntax: F = Function(@fun, any_data)
			%         F = Function(@fun)
			%
			%         F = Function(struct('p',2));
			%         F.setHandle(@(x)*F.Data.p)
			%
			% for more details, type "help Function"
			
			if nargin==0
				% nothing to do, an empty object will be automatically
				% constructed
			elseif nargin==1
				if isa(Handle, 'function_handle')
					F.Handle = Handle;
				else
					F.Data = Handle;
				end
			elseif nargin==2
				F.Handle = Handle;
				F.Data = Data;
			end
		end
		
		function out = feval(obj, x)
			% evaluates the function at a point
			
			% For performance reasons we really don't want to declare "feval"
			% as a method. Instead, remove this method, and rename the
			% "Handle" property to "feval". Then "obj.feval(x)" works
			% correctly without any overhead due to method call.
			out = obj.Handle(x);
		end
		
		function new = slice(obj, dims, values)
			% Slice a function through given coordinates

			narginchk(3, 3);
			
			% A very general implementation is to restore the vector "x" by
			% x(dims)=values, x(keep)=z, where "keep" are the indices which
			% are _not_ sliced.
			%
			% This implementation supports arbitrary nonlinear functions
			% specified as function handles, but also provides slicing of
			% OneNormFunction and InfNormFunction objects.
			new = Function(@(z) obj.Handle(obj.restore_sliced_x(z, dims, values)));
		end

	end
	
	methods (Hidden)
		
		function obj = setInternal(obj,name,value)
			%
			% overwrites internal property for the Function object
			% (internal function)
			%
			% If we want to add the internal property from outside of
			% this class (e.g. inside the PLCP solver) e.g.
			%
			% obj.Internal.name = value,
			%
			% use the syntax:
			%         obj.setInternal('name',value)
			%
			% DO NOT USE THIS METHOD UNLESS YOU PERFECTLY KNOW WHAT
			% YOU ARE DOING
			%
			
			narginchk(3, 3);
			
			if ~ischar(name)
				error('Name must be a string.');
			end
			
			obj.Internal.(name) = value;
			
		end
	end
	
	methods (Static, Hidden)
		function x = restore_sliced_x(z, dims, values)
			% Restores "x" by fixing x(dims)=values and x(other)=z with
			% "other = setdiff(1:nx, dims)
			
			nx = numel(z)+numel(dims);
			keep = setdiff(1:nx, dims);
			x = zeros(nx, 1);
			x(keep) = z(:);
			x(dims) = values;
		end
	end
	
end
