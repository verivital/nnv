function Q = uminus(P)
% Unitary minus. Q = -P.
% 
% @return Minus set
%

% allocate
Q(size(P)) = Polyhedron;

% revert the H- and V-representations, do not recompute vertices/halfspaces
% if they are not available in the source
for i=1:length(P)
    if P(i).hasHRep && ~P(i).hasVRep
        Q(i) = Polyhedron('H', [-P(i).A P(i).b], ...
            'He', [-P(i).Ae P(i).be]);
    elseif P(i).hasVRep && ~P(i).hasHRep
        Q(i) = Polyhedron('V', -P(i).V, 'R', -P(i).R);
    else
        Q(i) = Polyhedron('H',  [-P(i).A P(i).b], ...
            'He', [-P(i).Ae P(i).be], ...
            'V', -P(i).V, 'R',-P(i).R);
    end
end


end
