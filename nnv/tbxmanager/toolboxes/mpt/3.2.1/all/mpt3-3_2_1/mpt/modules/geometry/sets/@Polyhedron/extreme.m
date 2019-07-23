function sol = extreme(obj, y)
% EXTREME Compute extreme points of this polyhedron.
%
% -------------------------------------------------------------------
% sol = extreme(y) : Compute an extreme point in the direction y.
%
% Solve the optimization problem:
%
%   J(y) = max y'x s.t. x in Set
%
% and return an optimizer x.
%
% Parameters:
%  y   - Vector of length d
%
% Returns:
%  sol - Structure with fields
%   x           - Extreme point or [] if emptyset
%   exitflag    - mptopt return type
%   supp        - Support of the set in the direction x
%
%

global MPTOPTIONS
if isempty(MPTOPTIONS)
    MPTOPTIONS = mptopt;
end

narginchk(2, 2);

if numel(obj)>1
	% deal with arrays
	sol = obj.forEach(@(x) x.extreme(y));

else
	% compute the extreme point
    validate_realvector(y);
    y=y(:);
    if length(y) ~= obj.Dim,
        error('The vector x must be of length %i.', obj.Dim);
    end
    
    lp   = obj.optMat;
    lp.f = obj.buildCost(-y(:)).f;
	lp.quicklp = true;
    opt = mpt_solve(lp);
    sol.exitflag    = opt.exitflag;
    sol.x       = [];
    sol.supp    = [];
    switch sol.exitflag
        case MPTOPTIONS.UNBOUNDED
            sol.supp = inf;
        case MPTOPTIONS.OK
            sol.x       = opt.xopt(1:obj.Dim);
            sol.supp    = -opt.obj;
    end
end

end
