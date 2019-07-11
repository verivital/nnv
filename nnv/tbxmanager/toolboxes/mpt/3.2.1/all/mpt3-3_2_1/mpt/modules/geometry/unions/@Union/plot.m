function h = plot(varargin)
%
% Plotting of unions
%
% Plot a single union (or an array thereof):
%   U.plot()
% Plot a single union with options:
%   U.plot('opt1', val1, ...)
% Plot multiple unions:
%   plot(U1, U2)
% Plot multiple options with options:
%   plot(U1, 'opt1', val1, ..., U2, 'opt2', val2, ...)
% Return the handle:
%   h = plot(...)

% Implementation:
%
% Union/plot() is the main dispatcher which just parses input arguments and
% calls Union/plot_internal() for each detected union.
%
% Do not redefine the plot() method in derived classes. Instead, implement
% the plotting functionality in plot_internal(). Do not forget to declare
% this method as protected in the class constructor (see Union.m).

error(nargoutchk(0,1,nargout));

h = [];

% split input arguments into objects and corresponding options
[objects, options] = parsePlotOptions('Union', varargin{:});
if numel(objects)==0
	% no objects to plot
	return
end

prevHold = ishold;
if ~ishold, 
    newplot;
end
hold on

% plot each object
h = [];
idx = 1; % index such that each elements is plotted in different color
for i = 1:numel(objects)
	for j = 1:numel(objects{i})
		hj = plot_internal(objects{i}(j), idx, options{i}{:});
		idx = idx + objects{i}(j).Num;
		h = [h; hj];
	end
end

if ~prevHold,
    hold off;
end
if nargout==0
	clear h
end

end
