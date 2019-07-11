function T = triangulate(P)
%
% TRIANGULATE Create a triangulation.
%
% -------------------------------------------------------------------
% Description
% -------------------------------------------------------------------
%
%
%
% -------------------------------------------------------------------
% Syntax
% -------------------------------------------------------------------
%
% T = P.triangulate    Triangulate the polytope P
% T = triangulate(P)
%
% Notes:
%  - Can only triangulate bounded polytopes
%  - Uses qhull
%
% -------------------------------------------------------------------
%

% use P.forEach for arrays
error(P.rejectArray());

if ~P.isBounded || ~P.isFullDim || P.isEmptySet
    error('Only bounded, non-empty polyhedra in the full dimension can be triangulated.');
end
P.minVRep();
V = P.V;

% Triangulate V
%K = delaunayn(V);
K = delaunayn(V, {'Qt', 'Qbb', 'Qc', 'Qz'}); 
T(size(K,1))=Polyhedron;
for i=1:size(K,1)
    T(i) = Polyhedron(V(K(i,:),:));
end

% attach functions if defined
T.copyFunctionsFrom(P);

end
