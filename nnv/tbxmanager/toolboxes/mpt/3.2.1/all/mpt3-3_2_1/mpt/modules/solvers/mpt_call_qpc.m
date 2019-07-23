function R = mpt_call_qpc(S)
%
% header file to be inserted from XML source

global MPTOPTIONS
if isempty(MPTOPTIONS) && ~S.test
    MPTOPTIONS = mptopt;
end

if ~any(strcmpi(S.problem_type,{'LP','QP'}))
    error('mpt_call_qpc: QPC solver does not solve %s problems!',S.problem_type);
end

% for LP set H=0
if isempty(S.H)
    S.H = zeros(S.n);
end

% QPC does not accept sparse matrices
if issparse(S.H)
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

show_display = 0;

if any(strcmpi(S.solver,{'qpip','qpc.qpip'}))

    if S.test        
        % primal-dual predictor-corrector algorithm    
        [x, ef, L] = qpip(S.H,S.f,S.A,S.b,S.Ae,S.be,S.lb,S.ub,show_display);
    else
        mu = MPTOPTIONS.modules.solvers.qpip.mu;
        method = MPTOPTIONS.modules.solvers.qpip.method;
        [x, ef, L] = qpip(S.H,S.f,S.A,S.b,S.Ae,S.be,S.lb,S.ub,show_display,mu,method);
    end
else 
    % dual active-set algorithm (default)
    [x, ef, L] = qpas(S.H,S.f,S.A,S.b,S.Ae,S.be,S.lb,S.ub,show_display);
end

% test for NaN, or Inf
if any(isnan(x)) || any(isinf(x))
    x = zeros(S.n,1);
    ef = -1;
end

% must check feasibility
if ef==0 
    if S.m>1
        if S.test
            ve = any(S.A*x>S.b+1e-5);
        else
            ve = any(S.A*x>S.b+MPTOPTIONS.rel_tol);
        end
    else
        ve = false;
    end
    if S.me>1
        if S.test
            vq = ( norm(S.Ae*x-S.be,Inf) > 1e-5 );
        else
            vq = ( norm(S.Ae*x-S.be,Inf) > MPTOPTIONS.rel_tol );
        end
    else
        vq = false;
    end
    if any(ve) || vq
        % retry with interior point method
        if S.test
            [x, ef, L] = qpip(S.H,S.f,S.A,S.b,S.Ae,S.be,S.lb,S.ub,show_display);
        else
            mu = MPTOPTIONS.modules.solvers.qpip.mu;
            method = MPTOPTIONS.modules.solvers.qpip.method;
            [x, ef, L] = qpip(S.H,S.f,S.A,S.b,S.Ae,S.be,S.lb,S.ub,show_display,mu,method);
        end
        % check again feasibility
        if ef==0
            if S.m>1
                if S.test
                    ve = any(S.A*x>S.b+1e-5);
                else
                    ve = any(S.A*x>S.b+MPTOPTIONS.rel_tol);
                end
            else
                ve = false;
            end
            if S.me>1
                if S.test
                    vq = ( norm(S.Ae*x-S.be,Inf) > 1e-5 );
                else
                    vq = ( norm(S.Ae*x-S.be,Inf) > MPTOPTIONS.rel_tol );
                end
            else
                vq = false;
            end
            if any(ve) || vq                
                % infeasible
                ef = -1;
            end
        end
    end
end

R.xopt = x;
R.obj = 0.5*x'*S.H*x + S.f'*x;
R.lambda.ineqlin = L.inequality;
R.lambda.eqlin = L.equality;
R.lambda.lower = L.lowerbound;
R.lambda.upper = L.upperbound;

if ef==0
    R.how = 'ok';
    if S.test
        R.exitflag = 1;
    else
        R.exitflag = MPTOPTIONS.OK;
    end
elseif ef==-1
    R.how = 'infeasible';
    if S.test
        R.exitflag = 2;
    else
        R.exitflag = MPTOPTIONS.INFEASIBLE;
    end    
else
    R.how = 'unknown (possibly infeasible)';
    if S.test
        R.exitflag = 2;
    else
        R.exitflag = MPTOPTIONS.INFEASIBLE;
    end
    % since we don't know what the other flags mean, throw an error
    %error('mpt_call_qpc: other error code in "exitflag" from qpc.')
end
