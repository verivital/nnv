function sol = facetInteriorPoints(P)
%
% Compute points in the relative interior of each facet.
% 
% syntax: sol = facetInteriorPoints(P)
% 
% The output is a matrix composed of points listed rowwise where ith row
% corresponds to ith facet.
%

global MPTOPTIONS
if isempty(MPTOPTIONS)
    MPTOPTIONS = mptopt;
end

% use P.forEach() for arrays
error(P.rejectArray());

%% compute interior points via projection
% check if P is not empty
if P.isEmptySet
    sol = [];
    return;
end

% Can't compute points in the facets without the facets
P.minHRep();

% Compute interior point of the polyhedron
xp = P.interiorPoint;
% get the point
x = xp.x;

m = size(P.H,1);
sol = zeros(m, P.Dim);
for j=1:m
    n = P.H(j,1:end-1);
    
    % Project x onto each facet to see if it's interior
    alpha = (P.H(j,end) - n*x) / (n*n');
    y = x + alpha * n';
    
    % Test if y is in the interior of the facet
    s = P.H*[y;-1];
    s(j) = [];
    % test equalities as well
    t = P.He*[y;-1];
    if all(s < MPTOPTIONS.abs_tol) && norm(t,Inf)<MPTOPTIONS.rel_tol
        sol(j,:) = y';
    else
        % Compute an interior point on the facet by solving an LP
        s = P.interiorPoint(j);
        if isempty(s.x)
            error('Could not compute interior point in the %i''th facet.', j);
        end
        sol(j,:) = s.x';
    end
    
end

end
