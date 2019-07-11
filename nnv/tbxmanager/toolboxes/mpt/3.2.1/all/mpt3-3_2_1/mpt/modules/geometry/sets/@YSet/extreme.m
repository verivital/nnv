function sol = extreme(obj, x)
%
% ret = extreme(obj, x)
%
% Compute an extreme point of this set in the direction x.
%
% Parameters:
% x - Vector of length d
%
% Returns (as elements of ret)
% e    - An extreme point of C from the face , or NULL if empty.
% flag -  mptopt.INFEASIBLE
%         mptopt.UNBOUNDED
%         mptopt.OK
%

global MPTOPTIONS
if isempty(MPTOPTIONS)
    MPTOPTIONS = mptopt;
end

narginchk(2, 2);

% check vector x
validate_realvector(x);

% deal with arrays
if numel(obj)>1
	% return an array of structures
	sol = obj.forEach(@(elem) elem.extreme(x));
    return
end

% check dimension
if numel(x)~=obj.Dim
    error('The argument must have %i number of elements.', obj.Dim);
end

if any(size(x)~=size(obj.vars))
    x = transpose(x);
end

% if x was created out of the matrix, there might be symmetric terms, 
% the assign command checks for the compatibility of the vector
assign(obj.vars, x);

% make column vector out of any matrix
x = x(:);

model = obj.extr.model;
model.c = 0*model.c;
model.c(obj.extr.local) = -x;
s  = feval(model.solver.call,model);

% if we don't know if it is feasible or unbounded, retry with artificial bounds
if ismember(s.problem,[12, 15])
    model.lb = -MPTOPTIONS.infbound*ones(size(model.lb));
    model.ub = MPTOPTIONS.infbound*ones(size(model.ub));
    s = feval(model.solver.call,model);
    
    % if solution is feasible -> unbounded
    if s.problem == 0
        s.problem = 2;
    else
        % otherwise infeasible
        s.problem = 1;
    end

end

% get MPT flags
sol = yalmip2mptflag(s);

switch sol.exitflag
    case MPTOPTIONS.OK
        v = s.Primal;
        v = v(obj.extr.local(:));
        sol.x       = v;
        sol.supp    = x'*v;
        if sol.supp >= MPTOPTIONS.infbound
            sol.exitflag = MPTOPTIONS.UNBOUNDED;
            sol.supp = Inf;
        elseif sol.supp <= -MPTOPTIONS.infbound
            sol.exitflag = MPTOPTIONS.UNBOUNDED;
            sol.supp = -Inf;            
        end
    case MPTOPTIONS.INFEASIBLE
        sol.x = [];
        sol.supp = NaN;
    case MPTOPTIONS.UNBOUNDED;
        sol.x = s.Primal(obj.extr.local(:));
        sol.supp = Inf;
    otherwise
        error('Solver returned "%s" error when called from YALMIP.',sol.how);
end

end
