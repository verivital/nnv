function iMap = incidenceMap(P)
% INCIDENCEMAP Compute the incidence map of this polyhedron.
%
% iMap = incidenceMap()
%
% Map between which vertices are contained in which facets.
% 
% Note: This function computes both irredundant V and H-representations of
% the polyhedron and can be time consuming.
%
% Returns:
%  iMap  - Structure containing
%    V   - Vertices nV x d
%    R   - Rays     nR x d
%    H   - Facets   nF x d
%    He  - Affine hull 
%  incVH - nV x nF matrix. Element i,j is 1 if vertex i is contained in
%          facet j, 0 otherwise
%  incRH - nR x nF matrix. Element i,j is 1 if ray i is contained in
%          facet j, 0 otherwise
%
%

global MPTOPTIONS
if isempty(MPTOPTIONS)
    MPTOPTIONS = mptopt;
end

% deal with arrays
if numel(P)>1
	% return an array of structures
	iMap = P.forEach(@(elem) elem.incidenceMap);
    return
end

% pre-alocate output
iMap = struct('V', [], 'R', [], 'H', [], 'He', [], 'incRH', [], 'incVH', []);

% Need H-rep and V-rep
hRep = P.minHRep();
vRep = P.minVRep();

iMap.V  = vRep.V; 
iMap.R  = vRep.R;
iMap.H  = hRep.H; 
iMap.He = hRep.He;

nV = size(iMap.V,1);
nR = size(iMap.R,1);
nF = size(iMap.H,1);

% Vertex map is easy - just test for inclusion
iMap.incVH = sparse(abs(hRep.H*[vRep.V -ones(nV,1)]')' < MPTOPTIONS.abs_tol);

% A ray r is in a facet F if there exists a vertex v in F s.t. v+r in F 
iMap.incRH = spalloc(nR,nF,3*nR+3*nF);
for i=1:nF
  % Vertices in this facet
  v = vRep.V(iMap.incVH(:,i),:);
  
  % Test each ray
  for j=1:size(v,1)
    r = vRep.R + repmat(v(j,:),nR,1);
    iMap.incRH(:,i) = iMap.incRH(:,i) | abs(hRep.H(i,:)*[r -ones(nR,1)]')' < MPTOPTIONS.abs_tol;
  end
end

end
