function obj = computeVRep(obj)
% H-to-V conversion with possibly redundant V-rep output

global MPTOPTIONS

% deal with arrays
if numel(obj)>1
	if nargout==0
		obj.forEach(@computeVRep);
	else
		obj = obj.forEach(@computeVRep);
	end
	return
end

if obj.hasVRep
	% nothing to do
	return
elseif ~obj.hasHRep
	% empty set
	obj.hasVRep = true;
	return
elseif obj.isEmptySet()
	% empty set = empty vertices and rays
	obj.V_int = zeros(0, obj.Dim);
	obj.R_int = zeros(0, obj.Dim);
	obj.hasVRep = true;
	return
elseif obj.isFullSpace()
	% R^n = zero vertex and all basis vectors as rays
	obj.V_int = zeros(1, obj.Dim);
	obj.R_int = [eye(obj.Dim); -eye(obj.Dim)];
	obj.hasVRep = true;
	return
end

if false
    % vertex enumeration via duality
    %
    % see Borrelli, Bemporad, Morari: Predictive control for linear and
    % hybrid systems, Section 5.4.4, pp. 82-84 
    shifted = false;
    if ~obj.contains(zeros(obj.Dim, 1))
        % shift the polyhedron such that it contains the origin in its
        % interior
        if obj.isBounded()
            xint = obj.interiorPoint().x;
        else
            Rmax = 1; % bound the chebyradius for better numerical robustness
            sol = fast_chebyCenter(obj.H, obj.He, Rmax);
            xint = sol.x;
        end
        P = obj + (-xint);
        shifted = true;
    else
        P = obj;
    end
    D = P.dual();
    D.computeHRep();
    E = Polyhedron('H', D.H, 'He', D.He);
    F = E.dual();
    if shifted
        % shift back to the original location
        F = F + xint;
    end
    F.minVRep(); % kick out redundant vertices
    obj.V_int = F.V;
    obj.R_int = F.R;
    obj.hasVRep = true;
    return
end

done = false;
backupTried = false;

% work with minimal H-representations to improve numerics
obj.minHRep();

% shift the polytope such that it contains the origin in its interior --
% helps CDD a lot
xc = zeros(obj.Dim, 1);
if obj.isBounded()
    xc = obj.chebyCenter.x;
end
Ai = obj.A;
bi = obj.b - Ai*xc;
Ae = obj.Ae;
be = obj.be - Ae*xc;
A = [Ae; Ai]; % equalities must come first
B = [be; bi];

% Do vertex enumeration with CDD
while ~done
	try
		% change almost-zero elements to zero
		A(abs(A)<MPTOPTIONS.zero_tol) = 0;
		B(abs(B)<MPTOPTIONS.zero_tol) = 0;
		
		% round H-representation to certain number of decimal places
		roundfactor = 10^15;
		A = round(A*roundfactor) / roundfactor;
		B = round(B*roundfactor) / roundfactor;
		
		lin = 1:size(obj.He_int, 1);
		% call CDD
		try
			s = cddmex('extreme', struct('A',A,'B',B,'lin',lin));
		catch
			% if CDD fails, retry with Matlab version
			[s.V,s.R,s.A,s.B] = mpt_nlrs('extreme',obj);
        end
        % shift vertices back
		s.V = s.V + repmat(xc', size(s.V, 1), 1);
        
		if size(s.V,1) == 0 % This is a cone... we need an explicit vertex
			obj.V_int = [obj.V_int; zeros(1,obj.Dim)];
		else
			obj.V_int = [obj.V_int; s.V];
		end
		obj.R_int = [obj.R_int; s.R];
		done = true;
	catch
		if backupTried
			obj.V_int = [obj.V_int; zeros(1,obj.Dim)];
			done=true;
			% this error appears usually when plotting polyhedra,
			% it is therefore disabled to show at least something
			% error('Could not compute vertices : Numerical problems in CDD.');
		end
		backupTried = true;
		
		% Use trick from mpt2.6 and reduce the calculation precision
		% CDD sometimes fails to compute extreme points correctly
		% this happens often when slopes of two hyperplanes are too close
		% that's why we use a fixed-point arithmetics
		
		thull = Polyhedron('H', fix(obj.H_int * 1e5) / 1e5, ...
			'He', fix(obj.He_int*1e5) / 1e5);
		thull.minHRep();
		% Do vertex enumeration with CDD
		A = [thull.He_int(:,1:end-1); thull.H_int(:,1:end-1)];
		B = [thull.He_int(:,end); thull.H_int(:,end)];
	end
end

obj.hasVRep = true;

% unset obj.optMat since the H-representation might have changed
obj.optMat = [];

end
