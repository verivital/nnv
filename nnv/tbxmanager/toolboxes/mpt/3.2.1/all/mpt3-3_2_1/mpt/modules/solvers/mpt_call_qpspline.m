function R = mpt_call_qpspline(S)
%
% header file to be inserted from XML source

global MPTOPTIONS
if isempty(MPTOPTIONS) && ~S.test
    MPTOPTIONS = mptopt;
end

if ~strcmpi(S.problem_type,'QP')
    error('mpt_call_clp: QPspline solver does not solve %s problems!',S.problem_type);
end

% substitute
H = S.H;
f = S.f;
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

% store merged constraints
Am = A;
bm = b;


% remove equality constraints if any
Ae = S.Ae;
be = S.be;
me = S.me;
kept_rows.eq = 1:S.me;
if me>0
    % factorize Ae to get
    %  Ae(Br,Bc)*x(Bc) + Ae(Br,Nc)*x(Nc) = be(Br) % must be invertible mapping
    
    % check rank of Ae
    if S.test
        re = rank(Ae);
    else
        re = rank(Ae,MPTOPTIONS.rel_tol);
    end
    % check consistency
    if S.test
        rc = rank([Ae be]);
    else
        rc = rank([Ae be],MPTOPTIONS.rel_tol);
    end
    
    % for underdetermined system check linearly dependent rows
    if re<me || rc<me
        
        % if the right hand side is not linearly dependent, infeasible solution
        if rc>re
            R.xopt = zeros(S.n,1);
            R.obj = 0;
            R.lambda.ineqlin = [];
            R.lambda.eqlin = [];
            R.lambda.lower = [];
            R.lambda.upper = [];
            if S.test
                R.exitflag = 2;
            else
                R.exitflag = MPTOPTIONS.INFEASIBLE;
            end
            R.how = 'infeasible';
            return
        end
        
        while me ~= re
            % find linearly dependent rows
            [~,~,p] = lu(sparse([Ae be]),'vector');
            rd = p(re+1:end);
            
            % remove linearly dependent rows
            Ae(rd,:) = [];
            be(rd) = [];
            me = me-length(rd);
            kept_rows.eq(rd) = [];
            
            if S.test
                re = rank(full([Ae be]));
            else
                re = rank(full([Ae be]),MPTOPTIONS.rel_tol);
            end
        end
        
    end
    
    [Le,Ue,pe,qe] = lu(sparse(Ae),'vector');
    Br = pe(1:re); Bc = qe(1:re);
    Nr = pe(re+1:end); Nc = qe(re+1:end);
    % Nr must be empty -> otherwise there are depentent rows in S.Ae
    % which must be removed
    
    % substitute x(Bc) = C*x(Nc) + D
    Aebn = S.Ae(Br,Nc);
    %iAebb = inv(S.Ae(Br,Bc));
    beb = S.be(Br);
    
    % use factorized solution to compute C
    % C = -S.Ae(Br,Bc)\Aebn;
    Cl = -linsolve(full(Le),Aebn,struct('LT',true));
    C = linsolve(full(Ue(:,1:re)),Cl,struct('UT',true));
    
    % use factorized solution to compute D
    % D = S.Ae(Br,Bc)\beb;
    Dl = linsolve(full(Le),beb,struct('LT',true));
    D = linsolve(full(Ue(:,1:re)),Dl,struct('UT',true));
    
    Abc = A(:,Bc); Anc = A(:,Nc);
    
    % modify inequality constraints
    %A = -Abc*iAebb*Aebn + Anc;
    A = Abc*C + Anc;
    %b = b - Abc*iAebb*beb;
    b = b - Abc*D;
    
    % modify cost
    %H = S.H(Nc,Nc) + Aebn'*iAebb'*S.H(Bc,Bc)*iAebb*Aebn - Aebn'*iAebb'*S.H(Bc,Nc) - S.H(Nc,Bc)*iAebb*Aebn;
    H = S.H(Nc,Nc) + C'*S.H(Bc,Bc)*C + C'*S.H(Bc,Nc) + S.H(Nc,Bc)*C;
    %f = S.H(Nc,Bc)*iAebb*beb - Aebn'*iAebb'*S.f(Bc) + S.f(Nc) - Aebn'*iAebb'*S.H(Bc,Bc)*iAebb*beb;
    f = S.H(Nc,Bc)*D + C'*S.f(Bc) + S.f(Nc) + C'*S.H(Bc,Bc)*D;
    
end

% actual dimensions
[m,n] = size(A);

% transformation to form
% min 0.5 x'Mx - b'x
% st:  l< Ax <u

% artificial lower bound
if ~S.test
    l = -MPTOPTIONS.infbound*ones(m,1);
else
    l = -1e9*ones(m,1);
end

if S.test
    [R.xopt, lambda, R.obj, exitflag, R.how, iter, time] =  QPspline(H, f, l, A, b);
else
    [R.xopt, lambda, R.obj, exitflag, R.how, iter, time] =  QPspline(H, f, l, A, b, MPTOPTIONS.modules.solvers.qpspline);
end


% extract Lagrange multipliers
R.lambda.ineqlin = -lambda(1:S.m);
if ~isempty(S.lb)
    R.lambda.lower = zeros(S.n,1);
    R.lambda.lower(kept_rows.lb) = -lambda(S.m+1:S.m+numel(kept_rows.lb));
else
    R.lambda.lower = zeros(S.n,1);
end
if ~isempty(S.ub) && isempty(S.lb)
    R.lambda.upper = zeros(S.n,1);
    R.lambda.upper(kept_rows.ub) = -lambda(S.m+1:S.m+numel(kept_rows.ub));
elseif ~isempty(S.ub) && ~isempty(S.lb)
    R.lambda.upper = zeros(S.n,1);
    R.lambda.upper(kept_rows.ub) = -lambda(S.m+numel(kept_rows.lb)+1:S.m+numel(kept_rows.lb)+numel(kept_rows.ub));
else
    R.lambda.upper = zeros(S.n,1);
end


% if there were equalities, map back to original solution
if S.me>0
    xopt = zeros(S.n,1);
    xopt(Nc) = R.xopt;
    xopt(Bc) = C*R.xopt + D;
    R.xopt = xopt;
    % solve overdetermined system to get multipliers for equalities
    % H*x + f + Am'*lambda_ineq + Ae'*lambda_eq = 0
    lambda_eq = zeros(S.me,1);
    lambda_eq(kept_rows.eq) = -Ae'\(S.H*R.xopt + S.f - Am'*lambda);
    
    % extend multipliers
    R.lambda.eqlin = lambda_eq;
    R.obj = 0.5*R.xopt'*S.H*R.xopt + S.f'*R.xopt;
else
    R.lambda.eqlin = []; 
end


% correct multipliers for lower/upper bounds
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
R.lambda.lower(~activelb) = 0;
R.lambda.upper(~activeub) = 0;

% recast exitflag to MPT type
switch exitflag
    case 1
        if S.test
            R.exitflag = 1;
        else
            R.exitflag = MPTOPTIONS.OK;
        end
    case -1
        if S.test
            R.exitflag = 2;
        else
            R.exitflag = MPTOPTIONS.INFEASIBLE;
        end

    case {-2,-3,-4}
        if S.test
            R.exitflag = -1;
        else
            R.exitflag = MPTOPTIONS.ERROR;
        end        
end
    
    

end
