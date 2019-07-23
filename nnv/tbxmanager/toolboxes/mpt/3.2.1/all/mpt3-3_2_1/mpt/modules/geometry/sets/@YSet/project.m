function sol = project(obj, x)
%
% Compute the projection of the point x onto this set.
% Computes the closest point in this set to the point x:
%
% Parameters:
% x - Vector of length d
%
% Returns (as elements of ret)
% p    - Closest point to x contained in this set, or NULL if empty.
% flag -  1  if infeasible
%         2  if unbounded
%         12 if infeasible or unbounded (can't tell which)
%         0  if ok
%


global MPTOPTIONS
if isempty(MPTOPTIONS)
    MPTOPTIONS = mptopt;
end

narginchk(2, 2);

% check vector x
validate_realvector(x);

% reject arrays
error(obj.rejectArray());

% check dimension
if ~isequal(size(obj.vars), size(x))
	error('Input argument "x" must be a %dx%d matrix.', ...
		size(obj.vars, 1), size(obj.vars, 2));
end

% if x was created out of the matrix, there might be symmetric terms, 
% the assign command checks for the compatibility of the vector
assign(obj.vars, x);

% solve the problem via YALMIP
cost = norm(x - obj.vars,2)^2;
d = solvesdp(obj.constraints, cost, obj.opts);

% if we don't know if it is feasible or unbounded, retry with artificial bounds
if ismember(d.problem,[12, 15])
    F = obj.contraints + [ -MPTOPTIONS.infbound*ones(size(obj.vars)) <= obj.vars <= MPTOPTIONS.infbound*ones(size(obj.vars)) ];
    d = solvesdp(F, cost, obj.opts);
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
        sol.x = double(obj.vars);
        sol.dist = sqrt(double(cost));
    case MPTOPTIONS.INFEASIBLE,
        sol.x = [];
        sol.dist = NaN;
    case MPTOPTIONS.UNBOUNDED,
        sol.x = [];
        sol.dist = Inf;
    otherwise
        error('Solver returned "%s" error when called from YALMIP.',sol.how);
end
end
