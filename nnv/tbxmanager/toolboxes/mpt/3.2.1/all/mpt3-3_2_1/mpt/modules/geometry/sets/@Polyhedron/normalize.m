function obj = normalize(obj)
% NORMALIZE Return normalized polyhedron such that each facet i 
% a_i'*x<=b_i is scaled with ||a_i||_2 = 1.
% 
% P.normalize

global MPTOPTIONS
if isempty(MPTOPTIONS)
    MPTOPTIONS=mptopt;
end

% deal with arrays
if numel(obj)>1
    for i=1:numel(obj)
        obj(i).normalize;
    end
    return
end


if obj.hasHRep
    
    A = obj.A;
    b = obj.b;

    % 2-norm of each facet
    n = matNorm(A);

	% normalize 0'*x<=+/-b to 0'*x<= +/- sign(b)
	% (correct sign is important as not to mess with trivial infeasibility)
	ZeroRows = (n<MPTOPTIONS.zero_tol);
	n(ZeroRows) = 1;
	b(ZeroRows) = sign(b(ZeroRows));

	% normalize each half-space (0*x<=b will be left intact)
	nA = A ./ repmat(n,1,size(A,2));
	nb = b ./ n;
	obj.H_int=[nA, nb];

    % normalize also equalities if present
    if ~isempty(obj.He_int)
        Ae = obj.Ae;
        be = obj.be;
        
        % 2-norm of each equality
        ne = matNorm(Ae);
        
        % if polyhedron was properly constructed, n should not contain zero
        % values
        if any(abs(ne)<MPTOPTIONS.zero_tol)
            nAe=Ae;
            nbe=be;
        else
            nAe = Ae ./ repmat(ne,1,size(Ae,2));
            nbe = be ./ repmat(ne,1,size(be,2));
        end
        obj.He_int=[nAe, nbe];
    end

end

end
