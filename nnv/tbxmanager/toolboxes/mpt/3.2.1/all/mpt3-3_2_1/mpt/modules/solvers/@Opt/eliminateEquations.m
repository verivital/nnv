 function S = eliminateEquations(S)
%
% Remove the equations from the system
%
% If there are any equations on the parameter, then these are not removed
%
% Input format:
%   min 0.5 * x'*H*x + (pF*th + f)'*x + th'*Y*th + C*th + c
%   s.t.  A*x <= b + pB*th,  Ae*x = be + pE*th
%
% Output format:
%   min 0.5 * y'*H*y + (pF*th + f)'*y + th'*Y*th + C*th + c
%   s.t. A*y <= b + pB*th, 0 = be + pE*th
% x = recover.Y*y + recover.th*[th;1]
%
% Note that the optimizer will be the same, but the optimal value will change
%

global MPTOPTIONS

if isempty(MPTOPTIONS)
    MPTOPTIONS = mptopt;
end

% deal with arrays
if numel(S)>1
    for i=1:numel(S)
        S(i).eliminateEquations;
    end
    return
end

if ~strcmpi(S.problem_type,{'LP','QP','MILP','MIQP'})
    error('No equality removal for %s type of problems.',S.problem_type);
end

% check integer variables
if ~isempty(S.vartype)
    if any(S.vartype=='S') || any(S.vartype=='N')
        error('Semicontinuous (S) or semiinteger (N) variables are not supported in elimination procedure for equality constraints.');
    end
    
    % find continuous variables
    ind_b = find(S.vartype=='B');
    if isempty(ind_b)
        ind_b = zeros(0,1);
    end
    ind_i = find(S.vartype=='I');
    if isempty(ind_i)
        ind_i = zeros(0,1);
    end
    ind_c = find(S.vartype=='C');
    nr = numel(ind_c);

    if nr<=0
        error('Elimination cannot proceed because there are no continuous variables to eliminate.');
    end
else
    % index sets in case no vartype is specified
    ind_c = 1:S.n;
    ind_i = zeros(0,1);
    ind_b = zeros(0,1);
end

% for no equalities, return the same structure
if (S.me==0)
    disp('Problem does not contain equality constraints.');
    return;
end


% % remove -Inf/Inf from bounds
% if all(isinf(S.ub))
%     S.ub = [];
% else
%     S.ub(isinf(S.ub))=MPTOPTIONS.infbound;
% end
% if all(isinf(S.lb))
%     S.lb = [];
% else
%     S.lb(isinf(S.lb))=-MPTOPTIONS.infbound;
% end

% take merged inequality constraints
A = S.Internal.A;
b = S.Internal.b;
pB = S.Internal.pB;

% check rank of equality constraints
re = rank(full([S.Ae -S.pE S.be]));

% for underdetermined system check linearly dependent rows
me = S.me;
Ae = S.Ae;
be = S.be;
pE = S.pE;
kept_rows = 1:me;
if re<me
    while me ~= re
        % find linearly dependent rows Ae*x - pE*th = be
        [~,~,p,~] = lu(sparse([Ae -pE be]),'vector');
        rd = p(re+1:end);
        
        % remove linearly dependent rows
        Ae(rd,:) = [];
        be(rd) = [];
        if S.isParametric
            pE(rd,:) = [];
        end
        me = me-length(rd);
        kept_rows(rd) = [];
        
        re = rank(full([Ae -pE be]));
    end
end
S.Internal.removed_rows.eqlin = setdiff(1:me,kept_rows);

% extract part of Ae that corresponds to continuous/binary/integer variables
AeC = Ae(:,ind_c);
rec = rank(AeC,MPTOPTIONS.abs_tol);

% check rank of Ae
if rec==0
    error('Rank of equality constraint matrix "Ae" is equal zero, you cannot remove equalities.');
end

% factorize AeC to get an invertible mapping
% AeB*xB + AeI*xI + AeC(Br,Bc)*x(Bc) + AeC(Br,Nc)*x(Nc) = be(Br) + pE(Br,:)*th
[Le,Ue,pe,qe] = lu(sparse(AeC),'vector');
if rank(full(Ue(:,1:rec)),MPTOPTIONS.abs_tol)~=rec
    % if invertibility is not achieved, we need to factorize differently
    % try full factorization but we with different combination of
    % variables
    for i=1:S.n-rec
        [Le,Ue,pe] = lu(sparse(AeC(:,i:rec+i-1)),'vector');
        qe = 1:S.n;
        if rank(full(Ue(:,1:rec)),MPTOPTIONS.abs_tol)~=rec
            continue
        else
            break;
        end
    end
    if rank(full(Ue(:,1:rec)),MPTOPTIONS.abs_tol)~=rec
        error('EliminateEquations: Could not find invertible submatrix for removing equalities.');
    end
end
Br = pe(1:rec); Bc = qe(1:rec);
Nr = pe(rec+1:end); Nc = qe(rec+1:end);

% new index sets based on decomposition of continuous variables
ind_m = ind_c(Bc); % basic continuous variables
ind_n = ind_c(Nc); % non-basic continuous variables
if isempty(ind_n)
    ind_n=zeros(0,1);
end


% substitute x(Bc) = C1*xB + C2*yB + C3*x(Nc) + D1 + D2*th
% xB are binary variables, yB are new binary created from integer variables
AeB = Ae(Br,ind_b);
AeI = Ae(Br,ind_i);
Aebn = Ae(Br,ind_n);
%iAebb = inv(S.Ae(Br,ind_m));
beb = be(Br);
if S.isParametric
    pEb = pE(Br,:);
end

% use factorized solution to compute C1
if ~isempty(ind_i)
    T = S.Internal.T;
    t = S.Internal.t;
end
% C1 = -S.Ae(Br,ind_m)\AeB;
% C2 = -S.Ae(Br,ind_m)\[AeI*T];
Cl1 = -linsolve(full(Le(1:rec,:)),AeB,struct('LT',true));
C1 = linsolve(full(Ue(:,1:rec)),Cl1,struct('UT',true));
if ~isempty(ind_i)
    % count integers
    Cl2 = -linsolve(full(Le(1:rec,:)),AeI*T,struct('LT',true));
    C2 = linsolve(full(Ue(:,1:rec)),Cl2,struct('UT',true));
end

% use factorized solution to compute C3
% C3 = -S.Ae(Br,ind_m)\Aebn;
Cl3 = -linsolve(full(Le(1:rec,:)),Aebn,struct('LT',true));
C3 = linsolve(full(Ue(:,1:rec)),Cl3,struct('UT',true));

% use factorized solution to compute D1
% D1 = S.Ae(Br,ind_m)\beb;
if isempty(ind_i)
    % no integer variables
    Dl1 = linsolve(full(Le(1:rec,:)),beb,struct('LT',true));
    D1 = linsolve(full(Ue(:,1:rec)),Dl1,struct('UT',true));
else
    Dl1 = linsolve(full(Le(1:rec,:)),beb-AeI*t,struct('LT',true));
    D1 = linsolve(full(Ue(:,1:rec)),Dl1,struct('UT',true));    
end

if S.isParametric
    % use factorized solution to compute D2
    % D2 = S.Ae(Br,ind_m)\pE(Br,:);
    Dl2 = linsolve(full(Le(1:rec,:)),pEb,struct('LT',true));
    D2 = linsolve(full(Ue(:,1:rec)),Dl2,struct('UT',true));
end

% simplify substutition
AB = A(:,ind_b);
AI = A(:,ind_i);
Abc = A(:,ind_m); Anc = A(:,ind_n);

% modify inequality constraints
% alpha*[xB; yB] + beta*xN <= gamma + delta*theta
if isempty(ind_i)
    % no integers
    alpha = Abc*C1 + AB;
else
    alpha = [Abc*C1 + AB, Abc*C2 + AI*T];
end
beta = Abc*C3 + Anc;
S.A = [alpha beta];
if isempty(ind_i)
    % no integers
    S.b = b - Abc*D1;
else
    S.b = b - Abc*D1 - AI*t;
end
S.Internal.A = S.A; % update merged constraints stored internally
S.Internal.b = S.b; % update merged constraints stored internally
if S.isParametric
    S.pB = pB - Abc*D2;
    S.Internal.pB = S.pB; % update merged constraints stored internally
end

% sizes of variables
nxb = numel(ind_b);
if ~isempty(ind_i)
    nyb = size(T,2);
else
    nyb = 0;
end
nxn = numel(Nc);
S.vartype = [repmat('B',1,nxb+nyb),repmat('C',1,nxn)];

% modify cost for [xB; yB; xN]
if isempty(S.H)
    f = zeros(nxb+nyb+nxn,1);
    if isempty(ind_i)
        % no integers nyb=0
        f(1:nxb) = C1'*S.f(ind_m) + S.f(ind_b);
        f(nxb+1:nxb+nxn) = C3'*S.f(ind_m) + S.f(ind_n);

        c = S.f(ind_m)'*D1 + S.c;
        if S.isParametric
            pF = zeros(nxb+nxn,S.d);
            pF(1:nxb,:) = S.pF(ind_b,:) + C1'*S.pF(ind_m,:);
            pF(nxb+1:nxb+nxn,:) = C3'*S.pF(ind_m,:) + S.pF(ind_n,:);
            Y = S.pF(ind_m,:)'*D2 + S.Y;
            Cx = D1'*S.pF(ind_m,:) + S.f(ind_m)'*D2 + S.C;
        end
    else        
        f(1:nxb) = C1'*S.f(ind_m) + S.f(ind_b);            
        f(nxb+1:nxb+nyb) = T'*S.f(ind_i) + C2'*S.f(ind_m);
        f(nxb+nyb+1:nxb+nyb+nxn) = C3'*S.f(ind_m) + S.f(ind_n);            

        c = S.f(ind_m)'*D1 + S.c + S.f(ind_i)'*t;
        if S.isParametric
            pF = zeros(nxb+nyb+nxn,S.d);
            pF(1:nxb,:) = S.pF(ind_b,:) + C1'*S.pF(ind_m,:);
            pF(nxb+1:nxb+nyb,:) = T'*S.pF(ind_i,:) + C2'*S.pF(ind_m,:);
            pF(nxb+nyb+1:nxb+nyb+nxn,:) = C3'*S.pF(ind_m,:) + S.pF(ind_n,:);
            Y = S.pF(ind_m,:)'*D2 + S.Y;
            Cx = D1'*S.pF(ind_m,:) + S.f(ind_m)'*D2 + S.C + t'*S.pF(ind_i,:);
        end                
    end
else
    H = zeros(nxb+nyb+nxn);
    f = zeros(nxb+nyb+nxn,1);
    if isempty(ind_i)
        % no integers nyb=0
        H(1:nxb,1:nxb) = S.H(ind_b,ind_b) + C1'*S.H(ind_m,ind_m)*C1 + S.H(ind_b,ind_m)*C1 + C1'*S.H(ind_m,ind_b);
        H(1:nxb,nxb+1:nxb+nxn) = C1'*S.H(ind_m,ind_m)*C3 + S.H(ind_b,ind_m)*C3 + S.H(ind_b,ind_n) + C1'*S.H(ind_m,ind_n);
        H(nxb+1:nxb+nxn,1:nxb) = C3'*S.H(ind_m,ind_m)*C1 + C3'*S.H(ind_m,ind_b) + S.H(ind_n,ind_b) + S.H(ind_n,ind_m)*C1;
        H(nxb+1:nxb+nxn,nxb+1:nxb+nxn) = C3'*S.H(ind_m,ind_m)*C3 + S.H(ind_n,ind_n) + C3'*S.H(ind_m,ind_n) + S.H(ind_n,ind_m)*C3;
        
        f(1:nxb) = C1'*S.H(ind_m,ind_m)*D1 + 0.5*S.H(ind_m,ind_b)'*D1 + 0.5*S.H(ind_b,ind_m)*D1 + C1'*S.f(ind_m) + S.f(ind_b);
        f(nxb+1:nxb+nxn) = C3'*S.H(ind_m,ind_m)*D1 + 0.5*S.H(ind_m,ind_n)'*D1 + 0.5*S.H(ind_n,ind_m)*D1 + C3'*S.f(ind_m) + S.f(ind_n);

        c = 0.5*D1'*S.H(ind_m,ind_m)*D1 + S.f(ind_m)'*D1 + S.c;
        if S.isParametric
            pF = zeros(nxb+nxn,S.d);
            pF(1:nxb,:) = C1'*S.H(ind_m,ind_m)*D2 + 0.5*S.H(ind_m,ind_b)'*D2 + 0.5*S.H(ind_b,ind_m)*D2 + ...
                S.pF(ind_b,:) + C1'*S.pF(ind_m,:);
            pF(nxb+1:nxb+nxn,:) = C3'*S.H(ind_m,ind_m)*D2 + 0.5*S.H(ind_m,ind_n)'*D2 + 0.5*S.H(ind_n,ind_m)*D2 + ...
                C3'*S.pF(ind_m,:) + S.pF(ind_n,:);
            Y = 0.5*D2'*S.H(ind_m,ind_m)*D2 + S.pF(ind_m,:)'*D2 + S.Y;
            Cx = D1'*S.H(ind_m,ind_m)*D2 + D1'*S.pF(ind_m,:) + S.f(ind_m)'*D2 + S.C;
        end
    else
        H(1:nxb,1:nxb) = S.H(ind_b,ind_b) + C1'*S.H(ind_m,ind_m)*C1 + S.H(ind_b,ind_m)*C1 + C1'*S.H(ind_m,ind_b);
        H(1:nxb,nxb+1:nxb+nyb) = C1'*S.H(ind_m,ind_m)*C2 + S.H(ind_b,ind_i)*T + S.H(ind_b,ind_m)*C2 + C1'*S.H(ind_m,ind_i)*T;
        H(1:nxb,nxb+nyb+1:nxb+nyb+nxn) = C1'*S.H(ind_m,ind_m)*C3 + S.H(ind_b,ind_m)*C3 + S.H(ind_b,ind_n) + C1'*S.H(ind_m,ind_n);
        H(nxb+1:nxb+nyb,1:nxb) = C2'*S.H(ind_m,ind_m)*C1 + T'*S.H(ind_i,ind_b) + C2'*S.H(ind_m,ind_b) + T'*S.H(ind_i,ind_m)*C1;
        H(nxb+1:nxb+nyb,nxb+1:nxb+nyb) = T'*S.H(ind_i,ind_i)*T + C2'*S.H(ind_m,ind_m)*C2 + T'*S.H(ind_i,ind_m)*C2 + C2'*S.H(ind_m,ind_i)*T;
        H(nxb+1:nxb+nyb,nxb+nyb+1:nxb+nyb+nxn) = C2'*S.H(ind_m,ind_m)*C3 + T'*S.H(ind_i,ind_m)*C3 + T'*S.H(ind_i,ind_n) + C2'*S.H(ind_m,ind_n);
        H(nxb+nyb+1:nxb+nyb+nxn,1:nxb) = C3'*S.H(ind_m,ind_m)*C1 + C3'*S.H(ind_m,ind_b) + S.H(ind_n,ind_b) + S.H(ind_n,ind_m)*C1;
        H(nxb+nyb+1:nxb+nyb+nxn,nxb+1:nxb+nyb) = C3'*S.H(ind_m,ind_m)*C2 + C3'*S.H(ind_m,ind_i)*T + S.H(ind_n,ind_i)*T + S.H(ind_n,ind_m)*C2;
        H(nxb+nyb+1:nxb+nyb+nxn,nxb+nyb+1:nxb+nyb+nxn) = C3'*S.H(ind_m,ind_m)*C3 + S.H(ind_n,ind_n) + C3'*S.H(ind_m,ind_n) + S.H(ind_n,ind_m)*C3;
        
        f(1:nxb) = C1'*S.H(ind_m,ind_m)*D1 + 0.5*S.H(ind_m,ind_b)'*D1 + 0.5*S.H(ind_b,ind_m)*D1 + C1'*S.f(ind_m) + S.f(ind_b) + ...
            0.5*S.H(ind_i,ind_b)'*t + 0.5*S.H(ind_b,ind_i)*t + 0.5*C1'*S.H(ind_i,ind_m)'*t + 0.5*C1'*S.H(ind_m,ind_i)*t;
        f(nxb+1:nxb+nyb) = T'*S.H(ind_i,ind_i)*t + C2'*S.H(ind_m,ind_m)*D1 + 0.5*C2'*S.H(ind_i,ind_m)'*t + 0.5*C2'*S.H(ind_m,ind_i)*t + ...
            0.5*T'*S.H(ind_m,ind_i)'*D1 + 0.5*T'*S.H(ind_i,ind_m)*D1 + T'*S.f(ind_i) + C2'*S.f(ind_m);
        f(nxb+nyb+1:nxb+nyb+nxn) = C3'*S.H(ind_m,ind_m)*D1 + 0.5*S.H(ind_m,ind_n)'*D1 + 0.5*S.H(ind_n,ind_m)*D1 + C3'*S.f(ind_m) + S.f(ind_n) + ...
            0.5*C3'*S.H(ind_i,ind_m)'*t + 0.5*C3'*S.H(ind_m,ind_i)*t + 0.5*S.H(ind_i,ind_n)'*t + 0.5*S.H(ind_n,ind_i)*t;

        c = 0.5*D1'*S.H(ind_m,ind_m)*D1 + S.f(ind_m)'*D1 + S.c + ...
            0.5*t'*S.H(ind_i,ind_i)*t + 0.5*t'*S.H(ind_i,ind_m)*D1 + 0.5*D1'*S.H(ind_m,ind_i)*t + S.f(ind_i)'*t;
        if S.isParametric
            pF = zeros(nxb+nyb+nxn,S.d);
            pF(1:nxb,:) = C1'*S.H(ind_m,ind_m)*D2 + 0.5*S.H(ind_m,ind_b)'*D2 + 0.5*S.H(ind_b,ind_m)*D2 + ...
                S.pF(ind_b,:) + C1'*S.pF(ind_m,:);
            pF(nxb+1:nxb+nyb,:) = C2'*S.H(ind_m,ind_m)*D2 + 0.5*T'*S.H(ind_m,ind_i)'*D2 + 0.5*T'*S.H(ind_i,ind_m)*D2 + ...
                T'*S.pF(ind_i,:) + C2'*S.pF(ind_m,:);
            pF(nxb+nyb+1:nxb+nyb+nxn,:) = C3'*S.H(ind_m,ind_m)*D2 + 0.5*S.H(ind_m,ind_n)'*D2 + 0.5*S.H(ind_n,ind_m)*D2 + ...
                C3'*S.pF(ind_m,:) + S.pF(ind_n,:);
            Y = 0.5*D2'*S.H(ind_m,ind_m)*D2 + S.pF(ind_m,:)'*D2 + S.Y;
            Cx = D1'*S.H(ind_m,ind_m)*D2 + D1'*S.pF(ind_m,:) + S.f(ind_m)'*D2 + S.C + ...
                0.5*t'*S.H(ind_i,ind_m)*D1 + 0.5*t'*S.H(ind_m,ind_i)'*D1 + t'*S.pF(ind_i,:);
        end                
    end
    S.H = 0.5*(H+H');
end

S.f = f;
S.c = c;
if S.isParametric
    S.pF = pF;
    S.Y = Y;
    S.C = Cx;
end

% xM = C1*xB + C2*yB + C3*xN + D1 + D2*th
% xM = [C1 C2 C3]*[xB; yB; xN] + [D2 D1]*[th;1]

% Compute map from new variables back to old
% x = recover.Y*y + recover.th*[th;1]
%  [xB]   [ I  0  0] [xB]  [ 0   0] [th] 
%  [xI] = [ 0  T  0] [yB] +[ 0   t]*[1 ]
%  [xM]   [C1 C2 C3]*[xN]  [D2  D1]
%  [xN]   [ 0  0  I]       [ 0   0]

S.recover.Y = zeros(S.n,nxb+nyb+nxn);
S.recover.Y(ind_b,1:nxb) = eye(nxb);
if ~isempty(ind_i)
    S.recover.Y(ind_i,nxb+1:nxb+nyb) = T;
    S.recover.Y(ind_m,:) = [C1 C2 C3];
else
    S.recover.Y(ind_m,:) = [C1 C3];
end
S.recover.Y(ind_n,nxb+nyb+1:nxb+nyb+nxn) = eye(nxn);
S.recover.th = zeros(S.n,S.d+1);
if ~isempty(ind_i)
   S.recover.th(ind_i,S.d+1) = t;
end
if S.isParametric
    S.recover.th(ind_m,:) = [D2 D1];
else
    S.recover.th(ind_m,:) = D1;
end

% modify dimensions
S.n=nxb+nyb+nxn;
S.me=length(Nr);
S.m=size(S.A,1);

if S.me>0
    disp('EliminateEquations: The problem has been transformed, but some of the equality constraints could not have been removed.');
end

% remaining equality constraints on the parameter
% 0 = be + pE*th
% correct empty matrices with proper dimensions
S.Ae = zeros(S.me,S.n);
S.be = zeros(S.me, 1);
if S.isParametric
    S.pE = zeros(S.me,S.d);
end

% remaining equality constraints
if ~isempty(Nr)
    if isempty(ind_i)
        % no integers
        S.Ae = [Ae(Nr,ind_b)+Ae(Nr,ind_m)*C1, Ae(Nr,ind_m)*C3];
        S.be = be(Nr) - Ae(Nr,ind_m)*D1;
        if S.isParametric
            S.pE = pE(Nr,:);
        end
    else
        S.Ae = [Ae(Nr,ind_b)+Ae(Nr,ind_m)*C1, Ae(Nr,ind_i)*T+Ae(Nr,ind_m)*C2, Ae(Nr,ind_m)*C3];
        S.be = be(Nr) - Ae(Nr,ind_i)*t - Ae(Nr,ind_m)*D1;
        if S.isParametric
            S.pE = pE(Nr,:) - Ae(Nr,ind_m)*D2; 
        end
    end
end



% set up lb/ub
if ~isempty(S.lb)
    lb = zeros(S.n,1);
    lb(nxb+nyb+1:nxb+nyb+nxn) = S.lb(ind_n);
    S.lb = lb;
end
if ~isempty(S.ub)
    ub = ones(S.n,1);
    ub(nxb+nyb+1:nxb+nyb+nxn) = S.ub(ind_n);
    S.ub = ub;
end

end
