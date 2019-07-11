function sol = interiorPoint(P, facetI)
% INTERIORPOINT Compute an interior point of the polyhedron. Strictly
% interior if possible.
%
% ------------------------------------------------------------------
% sol = P.interiorPoint(facetI)
% 
% If facetI is non-empty, then computes a point int he interior of the
% facetI'th facet.
% Note: If interior points of all facets are desired, then calling
% facetInteriorPoint is faster.
% 
% Returns:
%  sol.x        - interior point, or [] if none exists
%  sol.isStrict - true if x is in the strict interior, false otherwise
%  sol.r        - if P has an Hrep, r is the radius of the inscribed ball, otherwise r=[]

global MPTOPTIONS
if isempty(MPTOPTIONS)
    MPTOPTIONS = mptopt;
end

if nargin < 2
    facetI = [];
end

%% allocate output
% must be here to always have the same ordering of fields
sol = struct('x', [], 'isStrict', [], 'r', []);

%% deal with arrays
if numel(P)>1
	% return an array of structures
	sol = P.forEach(@(elem) elem.interiorPoint(facetI));
    return
end

%% try to reuse a stored information
if isempty(facetI) && ~isempty(P.Internal.InternalPoint)
	sol = P.Internal.InternalPoint;
	return
end	

%% checks for polyhedra

% if facets provided, P must be irredundant
if ~isempty(facetI)
	if ~P.irredundantHRep
		error('Polyhedron must be in minimal representation when you want compute any interior point of its facets. Use "minHRep()" to compute the facets.');
	end
end

%% compute interior points

% try specific problems first for faster computation
if isempty(facetI)
    if P.hasVRep
        % Average vertices and rays
		if isempty(P.V) && ~isempty(P.R)
			% only rays
			sol.x = mean(normalize(P.R), 1)';
		else
			% vertices and maybe rays
			sol.x = mean(P.V,1)';
			if size(P.R,1) > 0
				sol.x = sol.x + norm(sol.x)*mean(normalize(P.R), 1)';
			end
		end
        
        V = P.V;
		if ~isempty(V)
			V = V(2:end,:) - repmat(V(1,:),size(V,1)-1,1);
		end
        R = P.R;
        
        sol.isStrict = false;
        if rank([V;R], MPTOPTIONS.abs_tol) == P.Dim
            sol.isStrict = true;
		end
		if isempty(V)
			% the set is unbounded
			sol.r = Inf;
		end
		P.Internal.InternalPoint = sol;
        return
    end
    
    % It's an affine set
    if size(P.He_int,1) > 0 && isempty(P.H_int)
        % only feasible sets
        if size(P.He,1)<=P.Dim
            sol.x = P.Ae\P.be;
            sol.isStrict = false;
            sol.r = inf;
			P.Internal.InternalPoint = sol;
			return
        end
    end
else
	% check if facet is vector of indices
	validate_indexset(facetI);
end

% must test for H-rep in case P is empty
if P.hasHRep
    % by default call chebyCenter for all problems
    chb = P.chebyCenter(facetI);
    if chb.exitflag == MPTOPTIONS.OK
        sol.x = chb.x;
        sol.r = chb.r;
        if 2*chb.r > MPTOPTIONS.region_tol
            sol.isStrict = true;
        else
            sol.isStrict = false;
        end
        if ~isempty(P.He_int)
            sol.isStrict = false;
		end
		if isempty(facetI)
			P.Internal.InternalPoint = sol;
		end
    end
end


end
