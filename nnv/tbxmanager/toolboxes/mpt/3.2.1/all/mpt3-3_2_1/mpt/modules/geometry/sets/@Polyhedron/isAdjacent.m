function [ts, iP, iQ] = isAdjacent(P,Q,fP,fQ)
%
% test if regions P and Q share a common facet - or if the facet to facet
% property holds
%
% inputs: regions P, Q
%         fP, fQ - (optional) indices of the facets to test
% outputs: ts - return status if P and Q are adjacent or not
%          iP - index of a facet of region P that is adjacent to iQ
%          iQ - index of a facet of region Q that is adjacent to iP
%

global MPTOPTIONS
if isempty(MPTOPTIONS)
   MPTOPTIONS = mptopt;
end

narginchk(2, 4);

if nargin < 3
    [ts, iP, iQ] = isNeighbor(P,Q);
elseif nargin < 4
    [ts, iP, iQ] = isNeighbor(P,Q,fP);
else
    [ts, iP, iQ] = isNeighbor(P,Q,fP,fQ);
end

% iP and iQ can have more elements which means that P or Q have
% equalities written as double-sided inequalities. In this case we take
% only the first index and that (although matlab seems to have no problem
% with checking multiple facets).
% see -> test_polyhedron_isadjacent_09_pass.m

if ts
    % check if facet of P and facet of Q are equal
    if ~isempty(iP)
        f1 = P.getFacet(iP(1));
    else 
        f1 = P;
    end
    if ~isempty(iQ)
        f2 = Q.getFacet(iQ(1));
    else
        f2 = Q;
    end
    if f1~=f2
        ts = false;
    end
end

end
 
