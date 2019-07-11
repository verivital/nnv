function R = mpt_call_gurobi(S)
%
% header file to be inserted from XML source

global MPTOPTIONS
if isempty(MPTOPTIONS) && ~S.test
    MPTOPTIONS = mptopt;
end

assert(~isequal(S.problem_type, 'LCP'), 'mpt_call_gurobi: GUROBI solver does not solve LCP problems!');

% A*x <= rhs
model.A = sparse([S.Ae; S.A]);
model.rhs = full([S.be; S.b]);
model.sense = char(['='*ones(S.me, 1); '<'*ones(S.m, 1)]); 

% lb <= x <= ub
if isempty(S.lb)
    model.lb = -Inf(S.n, 1);
else
    model.lb = S.lb;
end
if isempty(S.ub)
    model.ub = Inf(S.n, 1);
else
    model.ub = S.ub;
end

% types of variables (Binary, Integer, ...)
if isfield(S, 'vartype') && ~isempty(S.vartype)
    model.vtype = S.vartype;
end

% minimize x'*H*x + obj'*x
model.obj = S.f;
if isequal(S.problem_type(end-1:end), 'QP')
    % quadratic term for QP/MIQP problems
    model.Q = sparse(S.H*0.5);
end

% set options
if S.test
    opts.OutputFlag = 0; % no verbosity
else
    opts = MPTOPTIONS.modules.solvers.gurobi;
end
% prevent gurobi from returning the INF_OR_UNBD status:
% http://www.gurobi.com/documentation/5.6/reference-manual/infunbdinfo
opts.InfUnbdInfo = 1;

% solve
result = gurobi(model,opts);

switch result.status
    case 'OPTIMAL'
        R.how = 'ok';
        if S.test
            R.exitflag = 1;
        else
            R.exitflag = MPTOPTIONS.OK;
        end
        
        % objective value
        R.obj = result.objval;
        
        % primal optimizer
        R.xopt = result.x;
        
        % dual optimizer
        if ~isfield(result, 'pi')
            % no dual variables for mixed-integer problems
            lambda = NaN(size(model.A, 1), 1);
        else
            lambda = -result.pi;
        end
        
    case 'UNBOUNDED'
        R.how = 'unbounded';
        if S.test
            R.exitflag = 3;
        else
            R.exitflag = MPTOPTIONS.UNBOUNDED;
        end
        R.obj = -Inf; % unbounded minimization
        R.xopt = NaN(S.n, 1);
        lambda = NaN(size(model.A, 1), 1);
        
    case {'INFEASIBLE', 'INF_OR_UNBD'}
        R.how = 'infeasible';
        if S.test
            R.exitflag = 2;
        else
            R.exitflag = MPTOPTIONS.INFEASIBLE;
        end
        R.obj = Inf; % infeasible minimization
        R.xopt = NaN(S.n, 1);
        lambda = NaN(size(model.A, 1), 1);
        
    otherwise
        R.how = result.status;
        if S.test
            R.exitflag = -1;
        else
            R.exitflag = MPTOPTIONS.ERROR;
        end
        R.obj = Inf;
        R.xopt = NaN(S.n, 1);
        lambda = NaN(size(model.A, 1), 1);
        
end

% dual variables
R.lambda.ineqlin = lambda(S.me+1:S.me+S.m);
R.lambda.eqlin = lambda(1:S.me);
R.lambda.lower = NaN(S.n, 1);
R.lambda.upper = NaN(S.n, 1);

end
