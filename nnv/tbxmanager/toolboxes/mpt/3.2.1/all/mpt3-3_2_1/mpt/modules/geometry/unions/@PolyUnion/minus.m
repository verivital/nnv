function Un = minus(U,Q)
%
% Minkowski difference for PolyUnion object
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
%  The algorithm for efficiently computing the Minkowski difference between a
%  union of polytopes and a polytope is based on a talk by 
%  S. Rakovic and D. Mayne entitled "Constrained Control Computations"
%  It was the keynote address at the GBT Meeting, London, November, 2002.
%
%
% USAGE:
%   R=U-Q
%   R=minus(U,Q,Options)
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
% Un    - PolyUnion object updated with the Pontryagin difference
%
% see also PLUS, MLDIVIDE
%
%
% Copyright is with the following author(s):
% Revised 2012, Martin Herceg, Automatic Control Laboratory, ETH Zurich
%
% (C) 2005 Mario Vasak, FER, Zagreb
%          mario.vasak@fer.hr  
% (C) 2003 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch
% (C) 2003 Pascal Grieder, Automatic Control Laboratory, ETH Zurich,
%          grieder@control.ee.ethz.ch
% (C) 2003 Mato Baotic, Automatic Control Laboratory, ETH Zurich,
%          baotic@control.ee.ethz.ch

global MPTOPTIONS
if isempty(MPTOPTIONS)
    MPTOPTIONS = mptopt;
end

if ~isa(Q,'Polyhedron'),
    error('Second input argument must be a Polyhedron object');
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
if Q.isEmptySet || U.Num==0
    Un = PolyUnion(U.Set);
    Un.Internal = U.Internal;
    Un.Data = U.Data;
    return;
end

% MINKOWSKI DIFFERENCE ON UNIONS OF POLYHEDRA
%=============================================

% if we want to subtract full-dimensional polyhedron Q from any
% low-dimensional polyhedron- this will be empty, thus we remove any
% low-dim polyhedra first to simplify computations
if Q.isFullDim
    U.Set(~U.Set.isFullDim) = [];
end

% compute convex hull
Phull=U.convexHull;

% compute minkowski difference
PM = minus(Phull,Q);

% E is union of polytopes; sets inside Phull which are not covered by PA
E = Phull \ U.Set;

if numel(E)>0
    % discard low-dim polyhedra if all of U are full-dim
    E(~E.isFullDim)=[];

end
    
if numel(E)>0
    % flip Polytope
    QM = uminus(Q);
    
    % Minkowski addition on QM
    Emin = plus(E,QM);
    
    % Compute final polytopes
    R = PM \ Emin;
else
    R = PM;
end

% discard low-dim polyhedra if all of U are bounded
R(~R.isFullDim)=[];

% return new PolyUnion object
Un = PolyUnion(R);

% copy internal data including properties
Un.Internal = U.Internal;
Un.Data = U.Data;

% the new union does not have overlaps due to set-difference operation
if U.isFullDim
    Un.Internal.Overlaps = false;
end


end