function R = mpt_call_clp(S)
%
% header file to be inserted from XML source

global MPTOPTIONS
if isempty(MPTOPTIONS) && ~S.test
    MPTOPTIONS = mptopt;
end

if ~any(strcmpi(S.problem_type,{'LP','QP'}))
    error('mpt_call_clp: CLP solver does not solve %s problems!',S.problem_type);
end

% merge inequality constraints
A = S.A;
b = S.b;

% detect Inf boundaries
if S.test
    ilb = (S.lb==-Inf) | (S.lb<=-1e6);
    iub = (S.ub==Inf)  | (S.ub>=1e6);
else
    ilb = (S.lb==-Inf) | (S.lb<=-MPTOPTIONS.infbound);
    iub = (S.ub==Inf)  | (S.ub>=MPTOPTIONS.infbound);
end
% store kept rows
kept_rows.lb = find(~ilb);
kept_rows.ub = find(~iub);
if any(~ilb)
    % put ones at the positions where there is lb/ub
    Alb = zeros(nnz(~ilb),S.n);
    Alb(:,~ilb) = -eye(nnz(~ilb));
    A = [A; Alb];
    b = [b; -S.lb(~ilb)];
end
if any(~iub)
    Aub = zeros(nnz(~iub),S.n);
    Aub(:,~iub) = eye(nnz(~iub));
    A = [A; Aub];
    b = [b; S.ub(~iub)];
end

if ~S.test
    options=MPTOPTIONS.modules.solvers.clp;
    % function call by J. Loefberg
    [R.xopt, lambda, status] = clp(S.H, S.f, A, b, S.Ae, S.be, [], [], options);
else
    [R.xopt, lambda, status] = clp(S.H, S.f, A, b, S.Ae, S.be);
end


R.lambda.ineqlin = -lambda(S.me+1:S.me+S.m);
R.lambda.eqlin = -lambda(1:S.me);
if ~isempty(S.lb)
    R.lambda.lower = zeros(S.n,1);
    R.lambda.lower(kept_rows.lb) = -lambda(S.me+S.m+1:S.me+S.m+numel(kept_rows.lb));
else
    R.lambda.lower = zeros(S.n,1);
end
if ~isempty(S.ub) && isempty(S.lb)
    R.lambda.upper = zeros(S.n,1);
    R.lambda.upper(kept_rows.ub) = -lambda(S.me+S.m+1:S.me+S.m+numel(kept_rows.ub));
elseif ~isempty(S.ub) && ~isempty(S.lb)
    R.lambda.upper = zeros(S.n,1);
    R.lambda.upper(kept_rows.ub) = -lambda(S.me+S.m+numel(kept_rows.lb)+1:S.me+S.m+numel(kept_rows.lb)+numel(kept_rows.ub));
else
    R.lambda.upper = zeros(S.n,1);
end


if status==0,
    R.how = 'ok';
    if S.test
        R.exitflag = 1;
    else
        R.exitflag = MPTOPTIONS.OK;
    end
elseif status==1,
    R.how = 'infeasible';
    if S.test
        R.exitflag = 2;
    else
        R.exitflag = MPTOPTIONS.INFEASIBLE;
    end
elseif status==2
    R.how = 'unbounded';
    if S.test
        R.exitflag = 3;
    else
        R.exitflag = MPTOPTIONS.UNBOUNDED;
    end
else
    R.how = 'unknown error';
    if S.test
        R.exitflag = -1;
    else
        R.exitflag = MPTOPTIONS.ERROR;
    end
end

R.obj = S.f'*R.xopt;
if ~isempty(S.H)
    R.obj = R.obj + 0.5*R.xopt'*S.H*R.xopt;
end
