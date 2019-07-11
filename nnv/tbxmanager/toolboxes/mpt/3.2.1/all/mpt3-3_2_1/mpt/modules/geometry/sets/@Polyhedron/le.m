function tf = le(P,S)
% P less or equal than S

validate_polyhedron(S);

tf = S.contains(P);

end
