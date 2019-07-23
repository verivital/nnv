function Pnew = mtimes(P, S)
% MTIMES Multiply this polyhedron by a matrix, a scalar or a polyhedron.
%
% -------------------------------------------------------------------
% Description
% -------------------------------------------------------------------
%
% -------------------------------------------------------------------
% Syntax :
%    Pnew = mtimes(P, A)
%    Pnew = P*A
%    Pnew = A*P
%
% P is a polyhedron and A is a matrix or scalar
%
% Compute the affine map of a matrix A and this polyhedron.
%
% A*P = {A*x | x in P}
%
% Let A have n rows and d columns (the polyhedron must be in
% dimensional space).
%
% If n  < d then this is projection
% If n  > d then this is a lifting
% If n == d then this is rotation/skew
%
% If A is a scalar, then this function just scales the polyhedron.
%
% -------------------------------------------------------------------
% Syntax :
%    Pnew = mtimes(P, S)
%    Pnew = P*S
%    Pnew = S*P
%
% P and S are polyhedra
%
% Computes the cartesian product P x S
%

global MPTOPTIONS

% Determine the type of the arguments
if isa(S, 'Polyhedron')
    if numel(S)>1
        error('Only one polyhedron "S" can be provided.');
    end
    S_class = 'Polyhedron';
elseif isnumeric(S) && isscalar(S)
    validate_realvector(S);
    S_class = 'scalar';
elseif isnumeric(S)
    validate_realmatrix(S);
    S_class = 'matrix';
else
    error('This type of the object "%s" is not supported as argument for "S".', class(S));
end

if isa(P, 'Polyhedron')
    P_class = 'Polyhedron';
elseif isnumeric(P) && isscalar(P)
    validate_realvector(P);
    P_class = 'scalar';
elseif isnumeric(P)
    validate_realmatrix(P);
    P_class = 'matrix';
else
    error('This type of the object "%s" is not supported as argument for "P".', class(P));
end

% deal with arrays if P is polyhedron
if strcmpi(P_class,'Polyhedron') && numel(P)>1
    Pnew(size(P)) = Polyhedron;
    for i=1:numel(P)
        Pnew(i) = mtimes(P(i),S);
    end
    return;
end

switch [P_class S_class]
    case ['Polyhedron' 'Polyhedron']
        if S.isEmptySet,
            Pnew = Polyhedron(P);
            return;
        end
        if P.isEmptySet,
            Pnew = Polyhedron(S);
            return;
        end
        
        % Compute the cartesian product of P x S
        % For now - just compute the convex hull...
        P.minHRep(); S.minHRep();
		if isempty(P.He_int) && isempty(S.He_int)
			Pnew = Polyhedron(blkdiag(P.A, S.A), [P.b;S.b]);
		else
			Pnew = Polyhedron('H', [blkdiag(P.A, S.A) [P.b;S.b]], 'He', [blkdiag(P.Ae, S.Ae) [P.be;S.be]]);
		end
        %%% TODO : Compute polyhedral product without taking convex hull
        
    case {['Polyhedron' 'scalar'], ['scalar' 'Polyhedron']}
        if isnumeric(P),
            alpha = P;
            poly = S;
        else
            alpha = S;
            poly = P;
		end
		
		if abs(alpha) < MPTOPTIONS.abs_tol
			% scaling with zero produces a singleton
			%
			% we deal with this explicitly as to correctly support R^n
			% (depends on resolution of issue #93)

			% obey representation of the input
			if poly.hasHRep
				Pnew = Polyhedron([eye(poly.Dim); -eye(poly.Dim)], ...
					zeros(2*poly.Dim, 1));
			else
				Pnew = Polyhedron(zeros(1, poly.Dim));
			end
			
		elseif poly.hasHRep
			if isempty(poly.He_int)
				Pnew = Polyhedron(poly.A, alpha*poly.b);
			else
				Pnew = Polyhedron('H',[poly.A alpha*poly.b],'He', poly.He);
			end
        else
            Pnew = Polyhedron('V',alpha*poly.V, 'R', poly.R);
        end
    case ['Polyhedron' 'matrix']
        Pnew = P.invAffineMap(S);
    case ['matrix' 'Polyhedron']
        Pnew = S.affineMap(P);
end

end
