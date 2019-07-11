function R = mpt_call_qpoases(S)
%
% header file to be inserted from XML source

% requires QP to be formulated as:
%  min 0.5x'*H*x + f'*x
%   s.t.:  l <= x <= u
%          g1 <= G*x <= g2

global MPTOPTIONS
if isempty(MPTOPTIONS) && ~S.test
    MPTOPTIONS = mptopt;
end

if ~any(strcmpi(S.problem_type,{'QP','LP'}))
    error('mpt_call_qpoases: qpOASES solver does not solve %s problems!',S.problem_type);
end

% for LP we set H=0
if isempty(S.H)
    S.H = zeros(S.n);
end

% convert sparse matrices to full matrices
if issparse(S.H),
    S.H = full(S.H);
end
if issparse(S.f)
    S.f = full(S.f);
end
if issparse(S.A)
    S.A = full(S.A);
end
if issparse(S.b)
    S.b = full(S.b);
end
if issparse(S.Ae)
    S.Ae = full(S.Ae);
end
if issparse(S.be)
    S.be = full(S.be);
end

% equality constraints are treated as double-sided inequalities
G = [S.A; S.Ae];
g2 = [S.b; S.be];
% lower bound for the left hand side of inequalities g1 <= A*x<= g2 is
% set as very low number to achieve dimension match
if ~S.test
    g1 = [-MPTOPTIONS.infbound*ones(S.m,1); S.be];
else
    g1 = [-1e9*ones(S.m,1); S.be];
end
%The integer argument nWSR specifies the maximum number of working set
%recalculations to be performed during the initial homotopy (on output it contains the number
%of working set recalculations actually performed!)
nWSR = 5*(S.m+S.n); % this value is suggested by qpOASES manual


% call qpOASES
[R.obj, R.xopt, lambda, status]=qpOASES(S.H, S.f, G, S.lb,...
    S.ub, g1, g2, nWSR, S.x0);

% extract multipliers
if ~isempty(S.lb)
    if S.test
        activelb = (R.xopt < S.lb + 1e-4 );
    else
        activelb = (R.xopt < S.lb + MPTOPTIONS.rel_tol );
    end
else
    activelb = false(S.n,1);
end
if ~isempty(S.ub)
    if S.test
        activeub = (R.xopt > S.ub - 1e-4 );
    else
        activeub = (R.xopt > S.ub - MPTOPTIONS.rel_tol );
    end
else
    activeub = false(S.n,1);
end
R.lambda.lower = lambda(1:S.n);
R.lambda.lower(~activelb) = 0;
R.lambda.upper = -lambda(1:S.n);
R.lambda.upper(~activeub) = 0;
R.lambda.ineqlin = -lambda(S.n+1:S.n+S.m);
R.lambda.eqlin = -lambda(S.n+S.m+1:S.n+S.m+S.me);


switch status
    case 0
        R.how = 'ok';
        if S.test
            R.exitflag = 1;
        else
            R.exitflag = MPTOPTIONS.OK;
        end
    case 1
        R.how = ['Maximum number of working sets ', num2str(nWSR),' reached.'];
            %'You may change this value in S.options.nWSR.'];
        if S.test
            R.exitflag = -1;
        else
            R.exitflag = MPTOPTIONS.ERROR;
        end

    otherwise
        R.how = 'infeasible';
        if S.test
            R.exitflag = 2;
        else
            R.exitflag = MPTOPTIONS.INFEASIBLE;
        end

end
