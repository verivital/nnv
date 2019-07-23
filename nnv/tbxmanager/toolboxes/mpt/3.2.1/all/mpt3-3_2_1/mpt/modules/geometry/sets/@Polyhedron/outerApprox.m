function approx = outerApprox(obj)
%
% Computes the smallest axis-aligned hyperbox containing a given polyhedron

global MPTOPTIONS

if length(obj)>1
	for i = 1:length(obj)
		if nargout==1
			approx(i) = obj(i).outerApprox;
		else
			obj(i).outerApprox;
		end
	end
	return
end

d = obj.Dim;

if isfield(obj.Internal, 'lb') && isfield(obj.Internal, 'ub')
	% reuse stored information
	lb = obj.Internal.lb;
	ub = obj.Internal.ub;

elseif obj.hasVRep
	% Vrep is easy without rays, just take min/max of vertices
	if isempty(obj.R_int)
		lb = min(obj.V_int, [], 1)';
		ub = max(obj.V_int, [], 1)';
	else
		% resort to ConvexSet/outerApprox for unbounded polyhedra
		approx = outerApprox@ConvexSet(obj);
		lb = approx.Internal.lb;
		ub = approx.Internal.ub;
	end
	% update properties of the input object
	obj.Internal.lb = lb;
	obj.Internal.ub = ub;
	
elseif obj.isEmptySet()
    lb = Inf(obj.Dim, 1);
    ub = -Inf(obj.Dim, 1);
    obj.Internal.lb = lb;
    obj.Internal.ub = ub;
    
else
	% for Hrep we have to solve 2 LPs per each dimension
	H = obj.H_int;
	He = obj.He_int;
	f = zeros(1, d);
	LP.f = f;
	LP.A = H(:, 1:end-1);
	LP.b = H(:, end);
	LP.Ae = He(:, 1:end-1);
	LP.be = He(:, end);
	LP.lb = []; 
	LP.ub = [];
	LP.quicklp = true;
	
	lb = -Inf(d, 1);
	ub = Inf(d, 1);
	for i = 1:d
		% minimize
		LP.f = f;
		LP.f(i) = 1;
		sol = mpt_solve(LP);
		if sol.exitflag == MPTOPTIONS.OK
			lb(i) = sol.obj;
		end
		
		% maximize
		LP.f(i) = -1;
		sol = mpt_solve(LP);
		if sol.exitflag == MPTOPTIONS.OK
			ub(i) = -sol.obj;
		end
	end
	
	% update properties of the input object
	obj.Internal.lb = lb;
	obj.Internal.ub = ub;
end

% construct output arguments
if nargout==1
	% return the bounding hyperrectangle as a Polyhedron object
	approx = Polyhedron([eye(d); -eye(d)], [ub; -lb]);
	approx.irredundantHRep = true;
	approx.Internal.lb = lb;
	approx.Internal.ub = ub;
end
