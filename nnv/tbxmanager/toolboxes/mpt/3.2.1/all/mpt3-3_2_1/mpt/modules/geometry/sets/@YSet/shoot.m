function alpha = shoot(obj, x)
%
% Compute max alpha s.t. x*alpha in Set
%
% Parameters:
% x - Vector of length d
%
% Returns:
% alpha - Maximum value of alpha s.t. x*alpha is in set, or NaN if
% empty
%

global MPTOPTIONS
if isempty(MPTOPTIONS)
    MPTOPTIONS = mptopt;
end

narginchk(2, 2);

% check vector x
validate_realvector(x);

% deal with arrays
no = numel(obj);
if no>1
    alpha = NaN(size(obj));
    for i=1:no
        alpha(i) = obj(i).shoot(x);        
    end
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

% do not put obj.alpha directly into constraints because YALMIP
% changed its status to quadratic scalar which caused later an error when 
% detecting appropriate solver
a = obj.alpha;
% solve the problem via YALMIP
F = obj.constraints + [obj.vars(:) == x(:)*a];
d = solvesdp(F, -a, obj.opts);

% if we don't know if it is feasible or unbounded, retry with artificial bounds
if ismember(d.problem,[12, 15])
    F = F + [ -MPTOPTIONS.infbound*ones(size(obj.vars)) <= obj.vars <= MPTOPTIONS.infbound*ones(size(obj.vars)) ];
    F = F + [ -MPTOPTIONS.infbound <= a <= MPTOPTIONS.infbound ];
    d = solvesdp(F, -a, obj.opts);
    % if solution is feasible -> unbounded
    if d.problem == 0
        d.problem = 2;
    else
        % infeasible
        d.problem = 1;
    end
end

% get MPT flags
sol = yalmip2mptflag(d);

switch sol.exitflag
    case MPTOPTIONS.OK,
        alpha = double(a);
    case MPTOPTIONS.INFEASIBLE,
        alpha = NaN;
    case MPTOPTIONS.UNBOUNDED,
        alpha = Inf;
    otherwise
        error('Solver returned "%s" error when called from YALMIP.',sol.how);
end
end
