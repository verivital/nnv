function Pnew = homogenize(P, type)
%
% Compute a lifting of the polyhedron
%
% type = 'Hrep' or 'Vrep'
% Default = same as current rep
%

if nargin < 2 
    type = '';
end

%% deal with arrays
if length(P)>1
    Pnew(size(P)) = Polyhedron;
    for i=1:length(P)
        Pnew(i) = P(i).homogenize(type);
    end
    return;
end

if isempty(type)
    if P.hasHRep,
        type = 'hrep';
    else
        type = 'vrep';
    end
end


% empty polyhedron
if ~P.hasHRep && ~P.hasVRep
    Pnew = Polyhedron;
    return;
end


switch lower(type)
    
    case 'hrep'
        if ~P.hasHRep,
            P.minHRep();
        end
        Pnew = Polyhedron('H', [P.A -P.b zeros(size(P.H,1),1)], 'He', [P.Ae -P.be zeros(size(P.He,1),1)]);
    case 'vrep'
        if ~P.hasVRep,
            P.minVRep();
        end
        Pnew = Polyhedron('V', zeros(1,P.Dim+1), 'R', [P.V ones(size(P.V,1),1);P.R zeros(size(P.R,1),1)]);
    otherwise
        error('Unknown type %s', type);
end

end
