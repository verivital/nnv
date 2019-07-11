function [index, details] = locatePoint(U,x)
% Implementation of a graph search algorithm for a point location problem.
%
% The algorithm solves a point location problem for a convex union of
% non-overlapping polyhedra that comes from the PLCP solver. The
% implemented method is a graph traversal approach that uses the adjacency
% list returned from PLCP solver. If the explicit controller was generated
% using other parametric solver, the method is not applicable.
%
% Input: a point x in R^n
% Ouput: - index of a region where the region is located
%        - details with the number of operations needed
%

narginchk(2, 2);
% use U.forEach(@(u) u.contains(x)) to evaluate arrays
error(U.rejectArray());

% empty polyunion
if numel(U)==0 || U.Num<1
    index = [];
    details = [];
    return;
end

% check x
validate_realvector(x);
x = x(:);
if numel(x)~=U.Dim
    error('The point must have the same dimension as the union.');
end


% check if the adjacency list is present
if ~isfield(U.Internal,'adj_list')
    error('The union does not have an adjacency list. Please, use the polyunion output from PLCP solver that contains the adjacency list.');
end

% check the properties of the union
if isempty(U.Internal.Convex) || isempty(U.Internal.Overlaps) || ...
    isempty(U.Internal.Connected) || isempty(U.Internal.Bounded) || ...
    isempty(U.Internal.FullDim)        
    disp('Some of the required properties have not been determined for this union.')
    disp('The following properties are checked: convexity, overlaps, connectivity, boundedness, and full-dimensionality.');
    disp('This may take some time...');
end
if ~(U.isBounded && U.isFullDim && U.isConnected && U.isConvex && ~U.isOverlapping)
    error(['This method supports unions of polyhedra that are convex, non-overlapping, '... 
        'bounded, full-dimensional, connected, and come from PLCP solver with an adjacency list.']);
end

        
% call the graph traversal method
[index, details] = find_region(x,U.Set, U.Internal.adj_list);


end

function [index, details] = find_region(x,Parray,graph, index, details)
%
% For given point x, find a region index from the partition Parray of the
% explicit solution where the point lies.
%
% input:
%   x - value of parameter
%   Parray - an array of regions 
%   graph - represented by adjacency list
%

global MPTOPTIONS
if isempty(MPTOPTIONS)
    MPTOPTIONS=mptopt;
end

DEBUG=0;

if nargin<4
    index = 1;
    details.niter = 0;
    details.noper = 0;
    details.multiplications = 0;
    details.summations = 0;
    details.comparisons = 0;
end
if nargin<5
    details.niter = 0;
    details.noper = 0;
    details.multiplications = 0;
    details.summations = 0;
    details.comparisons = 0;
end

details.niter = details.niter + 1;

% extract polyhedron
P = Parray(index);

if DEBUG
    if index<=1
        plot(Parray); 
        hold on;
        text(x(1),x(2),'x');
    end
    xc = chebyCenter(P);
    text(xc.x(1),xc.x(2),num2str(index));
end  

% get direction
%P.normalize; % polytopes are normalized when returned from PLCP solver
d = P.A*x - P.b; % compute distance
direction = d > MPTOPTIONS.abs_tol;
pd = find(direction);

% count operations
n=P.Dim;
nv = size(P.H,1);
details.multiplications = details.multiplications + nv*n;
details.summations = details.summations + nv*n;
details.comparisons = details.comparisons + nv;
details.noper = details.multiplications + details.summations + details.comparisons;

if any(direction)
    % pick the direction with maximum distance over positive d
    [~,id]= max(d(pd));
    direction = false(size(direction));
    direction(pd(id))=true;
    
    % add comparisons
    details.comparisons = details.comparisons + numel(pd)-1;
    details.noper = details.multiplications + details.summations + details.comparisons;
    
    % next region to check
    index_new = graph{index}(direction);
    if length(index_new)>1
        % if there are more ways, pick the one which is not empty
        v = cellfun(@isempty,index_new);
        index_new = index_new(~v);
        if isempty(index_new)
            % point out of feasible area
            index = [];
            return;
        end
    end

    % if there are more choices, pick first    
    index_new=index_new{1};
    if numel(index_new)>1
        % if the facet neighbors to more regions, pick the first
        index_new=index_new(1);
    end
    
    if isempty(index_new)
        % point out of feasible area
        index = [];
    else
        % if x does not lie in this region, continue
        [index,details] = find_region(x,Parray,graph, index_new, details);
    end
    
end

end
