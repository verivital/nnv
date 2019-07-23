function Un = plus(U,Q)
%
% Minkowski summation for PolyUnion object
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
%  The algorithm for efficiently computing the Minkowski summation between a
%  union of polyhedra and a polytope 
%
% USAGE:
%   R=U+Q
%   R=plus(U,Q)
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% U   - PolyUnion object (union of polyhedra in the same dimension)
% Q   - Polytope 
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% U     - polytope (or polytope array) describing the Pontryagin difference
%

global MPTOPTIONS
if isempty(MPTOPTIONS)
    MPTOPTIONS = mptopt;
end

if ~isa(Q,'Polyhedron'),
    error('Second input argument must be a Polyhedron object!');
end
if numel(Q)~=1 
    error('Only single polyhedron "Q" is accepted.');
end
if numel(Q)>1
    if any(U.Dim~=[Q.Dim])
        error('The polyhedron array "Q" must be in the same dimension as the union.');
    end
end

% if Q is empty, create a new copy and quickly return
if Q.isEmptySet
    Un = PolyUnion(U.Set);
    Un.Internal = U.Internal;
    Un.Data = U.Data;
    return;
end


% MINKOWSKI SUMMATION ON UNIONS OF POLYHEDRA
%=============================================

% compute Minkowski sum for each polyhedron in the array
P = U.Set + Q;

% compute convex hull of the union
Phull=U.convexHull;

% compute minkowski summation for the convex hull
PM = plus(Phull,Q);

% compute the difference between the summation for the convex hull and P
T = PM \ P;

if numel(T)>0
    % discard low-dim polyhedra
    T(~T.isFullDim)=[];
        
    % R is final union
    R = PM \ T;
    
    % discard low-dim polyhedra
    R(~R.isFullDim)=[];
else
    R = PM;
end

% return new PolyUnion object
Un = PolyUnion(R);

% copy internal data including properties
Un.Internal = U.Internal;
Un.Data = U.Data;

% the new union does not have overlaps due to set difference operation
Un.Internal.Overlaps = false;

end