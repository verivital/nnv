function R = mpt_call_sedumi(S)
%
% header file to be inserted from XML source

global MPTOPTIONS
if isempty(MPTOPTIONS) && ~S.test
    MPTOPTIONS = mptopt;
end

if ~any(strcmpi(S.problem_type,{'LP','QP'}))
    error('mpt_call_clp: SEDUMI solver does not solve %s problems!',S.problem_type);
end

% merge inequality constraints
A=S.A;
b=S.b;

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
m = size(A,1);

if strcmpi(S.problem_type,'LP')
   
    % convert LP to form accepted by Sedumi
    %     min   c'*x
    % s.t.:   A*x = b
    %          x >= 0
   
    % Equality constraints Ae*x = be are extended in variables xp, xm
    % x = xp - xm
    % Inequality constraies A*x <= b are replaced by slack variable y to be
    % positive, -A*x-y = -b, y>=0
    
    % the final form from equalities
    % [Ae -Ae 0]*[xp; xm; y] = be
   
    % the final form from inequalities
    % [-A A -I]*[xp; xm; y] = -b
    %   xp, xm, y >= 0
       
    An = [-A A -eye(m)];
    bn = -b;
    if S.me>0
        An = [S.Ae -S.Ae zeros(S.me,m); An];
        bn = [S.be; bn];
    end
    
    % set dimension of nonnegativity constraints
    K.l = S.n+S.n+m;
    
    % objective function
    cn = [S.f; -S.f; zeros(m,1)];
    
    % call sedumi    
    if S.test
        % use default settings, with fid=0 (no verbosity)
        [z,dz,status] = sedumi(sparse(An),sparse(bn),sparse(cn),K,struct('fid',0));
    else
        % use global MPT settings
        opts = MPTOPTIONS.modules.solvers.sedumi;
        [z,dz,status] = sedumi(sparse(An),sparse(bn),sparse(cn),K,opts);
    end

    % recover original variables
    if ~isempty(z)
       R.xopt = full(z(1:S.n)-z(S.n+1:2*S.n));
    else
       R.xopt = zeros(S.n,1); 
    end
            
else
    % for QP we need to create second order cone
    %
    %     min t
    % s.t.   Ae*x = be
    %        A*x <= b
    %   x'*(H/2)*x + f'*x <= t
    %
    % Linear inequalities and equalities are rewritten in xp, xm, y
    % variables as in LP.
    % The quadratic constraint can be written as
    %
    %  ||      Qx            ||
    %  || 0.5*(1 + f'*x -t)  ||_2   <= 0.5*(1-f'*x+t)
    %  
    % which is passed to Sedumi via substitution
    %  u = [Q*x; 0.5*(1+f'*x-t)];
    %  v = [0.5*(1-f'*x+t)];
    % 
    % where H/2 = Q'*Q
    
    % factor H/2 = Q'*Q
    if S.test
        Q = cholinc(sparse(S.H/2),1e-8);
    else
        Q = cholinc(sparse(S.H/2),MPTOPTIONS.abs_tol);
    end
    % objective function in variables z=[t(1); xp(n); xm(n); y(m); v(1); u(n+1)]
    cn =[1; zeros(3*S.n+m+2,1)];
    
    % inequalities and equalities are extended in z variables as in LP case
    An = [zeros(m,1) -A A -eye(m) zeros(m,S.n+1) zeros(m,1)];
    bn = -b;
    if S.me>0
        An = [zeros(S.me,1) S.Ae -S.Ae zeros(S.me,m) zeros(S.me,S.n+1) zeros(S.me,1); An];
        bn = [S.be; bn];
    end
    
    % add second order cone constraints
    An = [An;
        zeros(S.n,1) Q -Q zeros(S.n,m) zeros(S.n,1) -eye(S.n,S.n+1);
       -0.5   S.f'/2  -S.f'/2  zeros(1,m) 0 [zeros(1,S.n) -1];
        0.5  -S.f'/2   S.f'/2  zeros(1,m) -1 zeros(1,S.n+1)];
    bn = [bn; zeros(S.n,1); -0.5; -0.5];
    
    % t, xp, xm, y >= 0
    % u, v belong to quadratic cone
    K.f = 1; % t is a free variable
    K.l = S.n+S.n+m; % nonnegative variables
    K.q = S.n+1+1; % v, u belong to quadratic cone v >= norm(u)
    
    % call sedumi in sparse format
    if S.test
        % use default settings, with fid=0 (no verbosity)
        [z,dz,status] = sedumi(sparse(An),sparse(bn),sparse(cn),K,struct('fid',0));
    else
        % use global MPT settings
        opts = MPTOPTIONS.modules.solvers.sedumi;
        [z,dz,status] = sedumi(sparse(An),sparse(bn),sparse(cn),K,opts);
    end
        
    % recover original variables in full format
    if ~isempty(z)
        R.xopt = full(z(2:S.n+1)-z(S.n+2:2*S.n+1));
    else
        R.xopt = zeros(S.n,1);
    end
        
end

% recover multipliers for both LP/QP case because the order of constraints
% remains the same (equalities, inequalities, lb, ub)
R.lambda.ineqlin = full(dz(S.me+1:S.me+S.m));
R.lambda.eqlin = -full(dz(1:S.me));
if ~isempty(S.lb)
    R.lambda.lower = zeros(S.n,1);
    R.lambda.lower(kept_rows.lb) = full(dz(S.me+S.m+1:S.me+S.m+numel(kept_rows.lb)));
else
    R.lambda.lower = zeros(S.n,1);
end
if ~isempty(S.ub) && isempty(S.lb)
    R.lambda.upper = zeros(S.n,1);
    R.lambda.upper(kept_rows.ub) = full(dz(S.me+S.m+1:S.me+S.m+numel(kept_rows.ub)));
elseif ~isempty(S.ub) && ~isempty(S.lb)
    R.lambda.upper = zeros(S.n,1);
    R.lambda.upper(kept_rows.ub) = full(dz(S.me+S.m+numel(kept_rows.lb)+1:S.me+S.m+numel(kept_rows.lb)+numel(kept_rows.ub)));
else
    R.lambda.upper = zeros(S.n,1);
end

% if the field is missing, we assume infeasibility
if ~isfield(status,'numerr')
    status.numerr = 3;
end
if ~isfield(status,'pinf')
    status.pinf = 0;
end
if ~isfield(status,'dinf')
    status.dinf = 0;
end


% check numerr first
switch status.numerr
    case 0
        R.how = 'ok';      
    case 1
        R.how = 'ok, but with numerical problems';
    case 2
        R.how = 'numerical problems';
    otherwise
        R.how = 'unknown problem';
end

% check primal, dual feasiblity
if status.pinf==status.dinf && any(status.numerr==[0,1])
    if S.test
        R.exitflag = 1;
    else
        R.exitflag = MPTOPTIONS.OK;
    end
elseif status.pinf==1,
    R.how = [R.how,', primal infeasible'];
    if S.test
        R.exitflag = 2;
    else
        R.exitflag = MPTOPTIONS.INFEASIBLE;
    end
elseif status.dinf==1
    R.how = [R.how,', dual infeasible'];
    if S.test
        R.exitflag = 2;
    else
        R.exitflag = MPTOPTIONS.INFEASIBLE;
    end
else
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
