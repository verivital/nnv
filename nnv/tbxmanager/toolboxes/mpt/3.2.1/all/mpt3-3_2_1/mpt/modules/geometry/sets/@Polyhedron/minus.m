function Pdiff = minus(P, S)
% MINUS Compute the minkowski difference of S with this polyhedron.
%
% -------------------------------------------------------------------------
% Case 1: x is a vector of length P.Dim
%
% Px = minus(P,x)
% Px = P.minus(x)
% Px = P - x
%
% Computes the Minkowski sum:
%
%   P-x = { y-x | y in P}
%
% Parameters:
%   x  - Vector of length P.Dim
%
% Returns:
%   Px - Minkowski sum of this polyhedron P and -x
%
% -------------------------------------------------------------------------
% Case 2: S is a polyhedron
%
% Pdiff = minus(P,S)
% Pdiff = P.minus(S)
% Pdiff = P - S
%
% Computes the Minkowski sum:
%
%  P-S= { x in P ~|~ x+w in P forall w in Q }
%
% Parameters:
%   S - Polyhedron of dimension P.Dim
%
% Returns:
%   Pdiff - Minkowski difference of P and S
%

global MPTOPTIONS
if isempty(MPTOPTIONS)
    MPTOPTIONS = mptopt;
end


% deal with arrays
if numel(P)>1
    Pdiff(size(P)) = Polyhedron;
    for i=1:numel(P)
        Pdiff(i) = minus(P(i),S);
    end
    return;
end

type = class(S);

switch type
    case 'Polyhedron'
        
        if numel(S)>1
            error('Only one polyhedron S is allowed.');
        end
                
        % if P is empty array
        if builtin('isempty',P)
            Pdiff = Polyhedron;
            return
        end
        
        if isEmptySet(S)
            Pdiff = P;
            return;
        end
        if P.Dim ~= S.Dim,
            error('P and S must have the same dimension'); 
        end
        
        % P-S = {x in P | x+w in P forall w in S}
        %%% TODO : Figure out how to compute Pontyargin diffs for V-rep
        P.minHRep();
        S.minHRep();
        
        % The affine hull of S must be a subset of the affine hull of P, or the
        % difference is empty.
        if rank([P.He;S.He], MPTOPTIONS.abs_tol) > rank(S.He, MPTOPTIONS.abs_tol)
            % empty polyhedron in the same dimension
			Pdiff = Polyhedron.emptySet(P.Dim);
            return
        end
        
        % special case P==S
        if P==S
            % empty polyhedron in the same dimension
            Pdiff = Polyhedron.emptySet(P.Dim);
            return;
        end
        
        % Subtract the support of S for each inequality of P
        Hn = [P.A P.b-S.support(P.A')];
        if any(Hn(:,end)==-Inf)
            % empty polyhedron in the same dimension
            Pdiff = Polyhedron.emptySet(P.Dim);
            return;
        end
        % remove remaining Inf rows
        Hn((Hn(:,end))==Inf,:) = [];
        Pdiff = Polyhedron('H', [P.H; Hn], 'He', P.He);
                       
    otherwise
        
        x = S(:);
        validate_realvector(x);
        if length(x) ~= P.Dim,
            error('Length of the vector must be P.Dim = %i.', P.Dim);
        end
        
        % Shift the polyhedron by x
        Pdiff = P.plus(-x);
end



end
