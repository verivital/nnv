function Q = uplus(P)
% Unitary plus. Q = +P.
% 
% @return the same polyhedron
%

% allocate
Q(size(P)) = Polyhedron;

% copy polyhedron
for i=1:length(P)
    Q(i) = Polyhedron(P(i));
end



end
