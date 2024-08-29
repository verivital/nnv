function result = verify_specification(reachSet, property)
% verify a based on the interection between the reach set and halfspaces defining the property
%   (assumed to be the un_robust_region or unsafe region to prove)
% 
% Syntax:
%    [result] = verify_speficication(reachSet, property)
%
% Inputs:
%   - reachSet: computed output set of neural network (e.g. 1x1 Star)
%   - property: cell array of cell array defining all conditions that outputSet must satisfy
%
% Output: (should we change this to sat, unsat or unknown?)
%   - result: 0 ->  property failed
%             1 ->  property satisfied
%             2 ->  unknown


    R = reachSet;
    nr = length(R);    % number of output sets (for approx should be 1)
    
    % Process property to verify
    if iscell(property) % created from vnnlib (one or multiple halfSpaces)
        property = property{1};
        property = property.Hg; % property transformed into a HalfSpace(s)
    end

    % Begin verification
    np = length(property);
    if np == 1 % only one halfspace
        for k = 1:nr
            Set = R(k);
            if ~isa(Set, "Star")
                Set = Set.toStar;
            end
            if isa(Set.V, 'gpuArray')
                Set = Set.changeDevice('cpu');
            end
            S = Set.intersectHalfSpace(property.G, property.g); % compute intersection with unsafe/not robust region
            if isempty(S)
                result = 1; % no intersection with unsafe region = safe (unsat)
            else 
                result = 2; % intersection with safe and unsafe region = unknown or unsafe
                break;
            end
        end
    else
        cp = 1; % current halfspace we are looking at
        result = 1; % start assuming property is unsat (no intersection)
        while cp <= np % multiple halfspaces, which means OR assertion
            for k = 1:nr % check every reach set vs OR property
                Set = R(k);
                if ~isa(Set, "Star")
                    Set = Set.toStar;
                end
                if isa(Set.V, 'gpuArray')
                    Set = Set.changeDevice('cpu');
                end
                S = Set.intersectHalfSpace(property(cp).G, property(cp).g); 
                if isempty(S)
                    continue; % does nothing, just need an statement, wanted to make this clear
                else
                    result = 2; %  unknown if approx, sat if exact
                    return;
                end
            end
            cp = cp+1;
        end
    end

end % close function

