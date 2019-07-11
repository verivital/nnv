function tf = isEmptySet(P)
%
% checks if given Polyhedron is empty
%
% Empty region is consired if it is not feasible and if the diameter of the
% chebyshev ball inscribed inside the polyhedron is less than "region_tol".
%

global MPTOPTIONS
if isempty(MPTOPTIONS)
    MPTOPTIONS = mptopt;
end

% Deal with arrays of polyhedra
tf = false(size(P));
for i=1:length(P)
    % compute interior point only if P(i).Internal.Empty property has not been set
	if P(i).Dim==0 || (isempty(P(i).H_int) && isempty(P(i).He_int) && ...
			isempty(P(i).V_int) && isempty(P(i).R_int))
		P(i).Internal.Empty = true;
	elseif isempty(P(i).Internal.Empty)
        ip = P(i).interiorPoint;
        if ~isempty(ip.r) && 2*ip.r<MPTOPTIONS.region_tol
            % note that low-dimensional polyhedra have radius=0
            % we need to check how far is the interior point from the
            % inequalities
            if ip.isStrict
                % full-dim
                P(i).Internal.Empty = true;
            else          
                % low-dim
                if P(i).hasVRep
                    % calculate the maximum distance to vertices
                    hn = P(i).V-repmat(x',2,1);
                    d = sqrt(sum(hn.*hn,2));
                    dmax = max(d);
                else
                    % calculate the maximum violation of inequalities
                    hn = P(i).H*[ip.x;-1];
                    dmax = max(hn);                    
                end
                if dmax>MPTOPTIONS.abs_tol
                    % maximum allowable tolerance is violated, region is empty
                    P(i).Internal.Empty = true;
                else
                    % no constraints violated, region is not empty
                    P(i).Internal.Empty = false;

                end
			end
		elseif any(isnan(ip.x))
            % NaN in interior point means empty polyhedron
            P(i).Internal.Empty = true;
		else
			P(i).Internal.Empty = isempty(ip.x);
        end
    end
    tf(i) = P(i).Internal.Empty;
end
end
