function obj = computeHRep(obj)
% V-to-H conversion with possibly redundant H-rep output

global MPTOPTIONS

% deal with arrays
if numel(obj)>1
	if nargout==0
		obj.forEach(@computeHRep);
	else
		obj = obj.forEach(@computeHRep);
	end
	return
end

if obj.hasHRep
	% nothing to do
	return
elseif ~obj.hasVRep
	% empty set
	obj.hasHRep = true;
	return
elseif obj.isFullSpace()
	% R^n
	Rn = Polyhedron.fullSpace(obj.Dim);
	obj.H_int = Rn.H;
	obj.He_int = Rn.He;
	obj.hasHRep = true;
	return
end

% compute Hrep

obj.minVRep();
if isempty(obj.R_int) && obj.Dim>1 && obj.isFullDim() && ...
		size(obj.V_int, 1)>=obj.Dim+1 && size(obj.V_int, 2)<=3
	% try convhulln first; requires following conditions to be met:
	% * no rays
	% * dimension at least 2
    % * the set is full-dimensional
	% * the set has at least d+1 vertices

	x0 = obj.interiorPoint().x;
	V = obj.V_int;
    K = convhulln(V); % {'QJ', 'Qx', 'Qs'}
	d = size(V, 2);
	s.A = zeros(size(K, 1), d);
	s.B = ones(size(K, 1), 1);
	s.lin = [];
    for i = 1:size(K, 1)
        % each row of K contains indices of vertices that lie on the i-th
        % facet
        P = V(K(i, :), :);
        % compute the normal vector and the offset of the facet
        W = [P, -ones(d, 1)];
        [AB, ~] = qr(W'); % qr() is much faster than null()
        a = AB(1:d, end);
        b = AB(end, end);
        
        % determine the sign
        if a'*x0>b
            a = -a;
            b = -b;
        end
        s.A(i, :) = a';
        s.B(i) = b;
    end

else
    % Do facet enumeration with CDD
    s = cddmex('hull', struct('V', obj.V_int, 'R', obj.R_int));
end

Hall  = [s.A s.B];
H = Hall; H(s.lin, :) = [];
He = Hall(s.lin, :);
obj.H_int = H;
obj.He_int = He;
obj.hasHRep = true;

% unset obj.optMat since the H-representation might have changed
obj.optMat = [];

end
