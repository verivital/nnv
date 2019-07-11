function sol = shoot(obj, r, x0)
% SHOOT Maximize along the ray alpha*r + x0 within this set.
%
% -------------------------------------------------------------------
% Description
% -------------------------------------------------------------------
%
% Solve the optimization problem:
%
%   J(x0,r) = max alpha s.t. alpha*r + x0 in Set
%
% and return the optimal value alpha.
%
% -------------------------------------------------------------------
% Syntax
% -------------------------------------------------------------------
%
% sol = shoot(r, x0)
%
% Parameters:
%  r        - Ray direction
%  x0       - [0] Ray origin
%
% Returns:
%  sol - Structure with fields
%   alpha   - Maximum length of ray, [] if emptyset, inf if unbounded
%   flag    - mptopt return type
%   flagTxt - mptopt return type text
%   x       - Extreme point of ray : alpha*r + x0
%

global MPTOPTIONS
if isempty(MPTOPTIONS)
    MPTOPTIONS = mptopt;
end

validate_realvector(r);
r=r(:);

% check dimension
D = [obj.Dim];
if any(D(1)~=D)
    error('The polyhedron array must be in the same dimension.');
end

if nargin<3
    x0 = zeros(obj(1).Dim,1);
else
    validate_realvector(x0);
end

% deal with arrays
if numel(obj)>1
	% return an array of structures
	sol = obj.forEach(@(elem) elem.shoot(r,x0));
    return
end

% prealocate output
sol = struct('exitflag', [], 'x', [], 'alpha', []);

% check sizes of vectors
if length(r) ~= obj.Dim,
    error('The ray vector "r"  must have the length of %i.', obj.Dim);
end
if length(x0) ~= obj.Dim,
    error('The vector "x0" must have the length of %i.', obj.Dim);
end


if obj.Dim>0
    % Re-arrange pre-computed set description matrices to:
    %  A*blkdiag(r,I) [alpha;lam] <= b - A*[x0;0]
    S.A  = [obj.optMat.A(:,1:obj.Dim)*r obj.optMat.A(:,obj.Dim+1:end)];
    S.b  = obj.optMat.b-obj.optMat.A(:,1:obj.Dim)*x0;
    S.Ae = [obj.optMat.Ae(:,1:obj.Dim)*r obj.optMat.Ae(:,obj.Dim+1:end)];
    S.be = obj.optMat.be-obj.optMat.Ae(:,1:obj.Dim)*x0;
    
    S.lb = [-inf;obj.optMat.lb(obj.Dim+1:end)];
    S.ub = [ inf;obj.optMat.ub(obj.Dim+1:end)];
    S.f  = [-1;zeros(size(S.Ae,2)-1,1)];
    
    res = mpt_solve(S);
        
else
    res.exitflag = MPTOPTIONS.ERROR;
end

sol.exitflag = res.exitflag;
sol.x = [];
sol.alpha = Inf;

switch sol.exitflag
    case MPTOPTIONS.OK,
        sol.alpha = -res.obj;
        sol.x     = sol.alpha*r + x0;
%         if sol.alpha>=MPTOPTIONS.infbound
%             sol.alpha = Inf;
%         end
%     case MPTOPTIONS.UNBOUNDED
%         sol.alpha = Inf;
%     otherwise
%         sol.alpha = [];
%         sol.x = [];
end

end
