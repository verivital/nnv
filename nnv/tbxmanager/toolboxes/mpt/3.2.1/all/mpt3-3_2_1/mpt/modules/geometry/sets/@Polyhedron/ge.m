function tf = ge(P,S)
% P greater or equal than S

validate_polyhedron(S);

tf = P.contains(S);

end
