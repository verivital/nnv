function h = fplot(obj, varargin)
%
% Union/fplot() plots functions defined over a union (or over a derived class)
%
% U.fplot() -- plots the only function
% U.fplot('fname') -- plots function of given name
% U.fplot('opt1', value1, ...) -- plots the only function with options
% U.fplot('fname', 'opt1', value1, ...) -- plots a given function with options
% h = U.fplot(...) -- returns handle of the plot
%
% Only allows plotting of a single union and a single function at each
% time.
%
% Options are the same as in ConvexSet/fplot

% Implementation notes:
%  * This is just a gateway function which delegates the actual plotting to
%    ConvexSet/fplot.
%  * Arrays of unions must be rejected (otherwise we would need to return
%    the handles as a cell array, which we don't do). Plot multiple
%    unions manually with "hold on".
%  * It handles PolyUnion objects as well (in fact, Union/fplot is very
%    general and supports arbitrary unions)
%  * Do NOT define fplot() methods in subclasses of Union!

error(nargoutchk(0,1,nargout));
error(obj.rejectArray());

if iscell(obj.Set)
	% plot each element of the set separately
	
	% hold the plot for the first element of the array
	prevHold = ishold;
	if ~ishold,
		newplot;
	end
	hold('on');

	h = [];
	dim = 0;
	for i = 1:numel(obj.Set)
		% tell ConvexSet/fplot to NOT set hold and axis, since doing so for
		% each element can be slow
		hi = obj.Set{i}.fplot(varargin{:}, 'array_index', i, ...
			'use_hold', false, 'use_3dview', false);
		h = [h; hi];
		dim = max(dim, obj.Set{i}.Dim);
	end
	
	% hold off at the end
	if ~prevHold
		hold('off');
	end
	if dim >= 2
		view(3);
		axis tight
	end

else
	% simpler syntax for arrays
	h = obj.Set.fplot(varargin{:});
end

if nargout==0
	clear h
end

end
