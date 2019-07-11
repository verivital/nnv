function R = mpt_call_lcp(S)
%
% header file to be inserted from XML source

global MPTOPTIONS

if isempty(MPTOPTIONS) && ~S.test
    MPTOPTIONS = mptopt;
end

if ~any(strcmpi(S.problem_type,{'LP','QP','LCP'}))
    error('mpt_call_lcp: LCP solver does not solve %s problems!',S.problem_type);
end

% convert sparse matrices to full matrices
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

%% for LP, QP we need to transform to LCP form
if any(strcmpi(S.problem_type,{'LP','QP'}))
            
    % if no H is present
    if isempty(S.H)
        S.H = zeros(S.n);
    end
    H = S.H;
    f = S.f;
    
    % merge inequality constraints
    A=S.A;
    b=S.b;
    % detect Inf boundaries
	if ~S.test && (~isempty(S.lb) || ~isempty(S.ub))
        ilb = (S.lb==-Inf) | (S.lb<=-MPTOPTIONS.infbound);
        iub = (S.ub==Inf)  | (S.ub>=MPTOPTIONS.infbound);
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
	else
		kept_rows.lb = [];
		kept_rows.ub = [];
	end
    % store merged constraints
    Am = A;
    bm = b;
    
    % remove equality constraints if any
    kept_rows.eq = 1:S.me;
    Ae = S.Ae;
    be = S.be;
    me = S.me;
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

		if re>=S.n
			% overdetermined system, no degrees of freedom
			R.xopt = S.Ae\S.be;
			
			% check constraints
			eqc = ( norm(S.Ae*R.xopt-S.be,Inf) <= MPTOPTIONS.rel_tol );
			if ~isempty(S.b)
				ineqc = ( S.A*R.xopt <= S.b + MPTOPTIONS.rel_tol );
			else
				ineqc = true;
			end
			
			% multipliers
			% H*x + f + A'*lambda_ineq + Ae'*lambda_eq = 0
			R.lambda.ineqlin = zeros(S.m,1);
			R.lambda.eqlin = -transpose(S.Ae)\(H*R.xopt + S.f);
			R.lambda.lower = zeros(S.n,1);
			R.lambda.upper = zeros(S.n,1);
			
			if eqc && all(ineqc)
				% feasible solution
				if S.test
					R.exitflag = 1;
				else
					R.exitflag = MPTOPTIONS.OK;
				end
				R.how = 'ok';
			else
				% infeasible
				if S.test
					R.exitflag = 2;
				else
					R.exitflag = MPTOPTIONS.INFEASIBLE;
				end
				R.how = 'infeasible';
			end
			R.obj = 0.5*R.xopt'*H*R.xopt + S.f'*R.xopt + S.c;
			return
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
        
        if ~isempty(Ae)
            
            [Le,Ue,pe,qe] = lu(sparse(Ae),'vector');
            
            if S.test
                rk = rank(full(Ue(:,1:re)));
            else
                rk = rank(full(Ue(:,1:re)),MPTOPTIONS.abs_tol);
            end
            if rk~=re
                % if invertibility is not achieved try reduced echelon
                % elimination
                
                % find linear independent columns of Ae
                if S.test
                    [~,jb]=rref(Ae,1e-4);
                else
                    [~,jb]=rref(Ae,MPTOPTIONS.abs_tol);
                end
                [Le,Ue,pe] = lu(sparse(Ae(:,jb)),'vector');
                % update column selection
                qe = [jb, setdiff(1:S.n,jb)];
                
                % here check rank with default tolerance ->
                % test_opt_eliminateEquations_12_pass
                rk = rank(full(Ue(:,1:re)));
                if rk~=re
                    % all possibilities failed, don't know what to do more, for now just throw error
                    error('mpt_call_lcp: Could not find invertible submatrix for removing equalities.');
                end
                
            end
            Br = pe(1:re); Bc = qe(1:re);
            Nr = pe(re+1:end); Nc = qe(re+1:end);
            % Nr must be empty -> otherwise there are depentent rows in Ae
            % which must be removed
            
            % substitute x(Bc) = C*x(Nc) + D
            Aebn = Ae(Br,Nc);
            %iAebb = inv(Ae(Br,Bc));
            beb = be(Br);
            
            % use factorized solution to compute C
            % C = -Ae(Br,Bc)\Aebn;
            Cl = -linsolve(full(Le),Aebn,struct('LT',true));
            C = linsolve(full(Ue(:,1:re)),Cl,struct('UT',true));
            
            % use factorized solution to compute D
            % D = Ae(Br,Bc)\beb;
            Dl = linsolve(full(Le),beb,struct('LT',true));
            D = linsolve(full(Ue(:,1:re)),Dl,struct('UT',true));
            
            Abc = A(:,Bc); Anc = A(:,Nc);
            
            % modify inequality constraints
            %A = -Abc*iAebb*Aebn + Anc;
            A = Abc*C + Anc;
            %b = b - Abc*iAebb*beb;
            if ~isempty(b)
                b = b - Abc*D;
            else
                b = [];
            end
            
            % modify cost
            %H = S.H(Nc,Nc) + Aebn'*iAebb'*S.H(Bc,Bc)*iAebb*Aebn - Aebn'*iAebb'*S.H(Bc,Nc) - S.H(Nc,Bc)*iAebb*Aebn;
            H = S.H(Nc,Nc) + C'*S.H(Bc,Bc)*C + C'*S.H(Bc,Nc) + S.H(Nc,Bc)*C;
            %f = S.H(Nc,Bc)*iAebb*beb - Aebn'*iAebb'*S.f(Bc) + S.f(Nc) - Aebn'*iAebb'*S.H(Bc,Bc)*iAebb*beb;
            f = S.H(Nc,Bc)*D + C'*S.f(Bc) + S.f(Nc) + C'*S.H(Bc,Bc)*D;
            
        else
           % matrix Ae has been completely eliminated 
           S.me = 0;
           S.Ae = zeros(0,S.n);
           S.be = [];
        end
    end
    
    % actual dimensions
    if S.test
        r = rank(A);
    else
        r = rank(A,MPTOPTIONS.abs_tol);
    end
    %m = size(A,1);
    n = size(A,2);
    if r<n
        % express x as a difference of two positive numbers, i.e. x = x+ - x-
        
        % new objective function
        Hnew = [H -H; -H H];
        fnew = [f; -f];
        
        % new constraints Anew*y >= bnew
        Anew = [-A A];
        bnew = -b;
        
        % catch error when allocating on extra large dimensions
        try
            % create M,q for LCP
            M = [Hnew -Anew'; Anew zeros(size(Anew,1))];
            q = [fnew; -bnew];
            
            % solve LCP
            % dimension of LCP to solve is "2*size(H,1)+size(A,1)"
            if S.test
                [z,w,basis,exfl] = lcp(M, q);
            else
                if ~isempty(S.routine)
                    lcpopt = MPTOPTIONS.modules.solvers.lcp;
                    lcpopt.routine = S.routine;
                    [z,w,basis,exfl] = lcp(M, q, lcpopt );
                else
                    [z,w,basis,exfl] = lcp(M, q, MPTOPTIONS.modules.solvers.lcp );
                end
            end
        catch M1
            % display what happened wrong
            disp(M1.message);            
            z = zeros(size(H,1)+size(Anew,1),1);
            w = z;
            % return error status 
            exfl = -4;                      
        end
        
        % recover solution
        R.xopt = z(1:n) - z(n+1:2*n);

        % recover mulpliers
        lambda_ineq = z(2*n+1:end);
        R.lambda.ineqlin = lambda_ineq(1:S.m);
        if ~isempty(S.lb)
            R.lambda.lower = zeros(S.n,1);
            R.lambda.lower(kept_rows.lb) = lambda_ineq(S.m+1:S.m+numel(kept_rows.lb));
        else
            R.lambda.lower = zeros(S.n,1);
        end
        if ~isempty(S.ub) && isempty(S.lb)
            R.lambda.upper = zeros(S.n,1);
            R.lambda.upper(kept_rows.ub) = lambda_ineq(S.m+1:S.m+numel(kept_rows.ub));
        elseif ~isempty(S.ub) && ~isempty(S.lb)
            R.lambda.upper = zeros(S.n,1);
            R.lambda.upper(kept_rows.ub) = lambda_ineq(S.m+numel(kept_rows.lb)+1:S.m+numel(kept_rows.lb)+numel(kept_rows.ub));
        else
            R.lambda.upper = zeros(S.n,1);
        end
        
        
    else        
        % factorize A to get
        %  -A(B,:)*x = -b(B) + y  % must be invertible mapping
        %  -A(N,:)*x >= -b(N)
        %          y >= 0
        [L,U,p] = lu(sparse(-A),'vector');
        B = p(1:n);
        N = p(n+1:end);
        
        % substitute
        % use factorized solution to compute inv(A(B,:))
        iAbl = -linsolve(full(L(1:n,:)),eye(n),struct('LT',true));
        iAb = linsolve(full(U),iAbl,struct('UT',true));
        %iAb = inv(A(B,:)); 
        bb = b(B);
        An = A(N,:);
        bn = b(N);
        
        % form new objective function
        Hnew = iAb'*H*iAb;
        fnew = -iAb'*H*iAb*bb - iAb'*f;
        
        % new constraints Anew* y>= bnew
        Anew = An*iAb;
        if isempty(bn)
            bnew = [];
        else
            bnew = Anew*bb - bn;
        end
        
        % catch error when allocating on extra large dimensions
        try
            % create M,q for LCP
            M = [Hnew -Anew'; Anew zeros(size(Anew,1))];
            q = [fnew; -bnew];
                    
            % solve LCP
            % dimension of LCP to solve is "size(H,1)+size(Anew,1)"
            if S.test
                [z,w,basis,exfl] = lcp(M, q);
            else
                if ~isempty(S.routine)
                    lcpopt = MPTOPTIONS.modules.solvers.lcp;
                    lcpopt.routine = S.routine;
                    [z,w,basis,exfl] = lcp(M, q, lcpopt );
                else
                    [z,w,basis,exfl] = lcp(M, q, MPTOPTIONS.modules.solvers.lcp );
                end
            end
        catch M2
            % display what happened wrong
            disp(M2.message);
            z = zeros(size(H,1)+size(Anew,1),1);
            w = z;
            % return error status 
            exfl = -4;
        end
        
        % recover variables
        xopt = z(1:n);
        lambda_ineq = zeros(length(bm),1);
        lambda_ineq(B) = w(1:n);
        lambda_ineq(N) = z(n+1:end);        
                
        % recover solution
        R.xopt = iAb*(bb - xopt);
        % recover mulpliers
        R.lambda.ineqlin = lambda_ineq(1:S.m);
        if ~isempty(S.lb)
            R.lambda.lower = zeros(S.n,1);
            R.lambda.lower(kept_rows.lb) = lambda_ineq(S.m+1:S.m+numel(kept_rows.lb));
        else
            R.lambda.lower = zeros(S.n,1);
        end
        if ~isempty(S.ub) && isempty(S.lb)
            R.lambda.upper = zeros(S.n,1);
            R.lambda.upper(kept_rows.ub) = lambda_ineq(S.m+1:S.m+numel(kept_rows.ub));
        elseif ~isempty(S.ub) && ~isempty(S.lb)
            R.lambda.upper = zeros(S.n,1);
            R.lambda.upper(kept_rows.ub) = lambda_ineq(S.m+numel(kept_rows.lb)+1:S.m+numel(kept_rows.lb)+numel(kept_rows.ub));
        else
            R.lambda.upper = zeros(S.n,1);
        end

        
    end
    
    
    % if equalities were present, map back to original variables
    if S.me>0
        xopt = zeros(S.n,1);
        xopt(Nc) = R.xopt;
        %xopt(Bc) = iAebb*(-Aebn*R.xopt + beb);
        xopt(Bc) = C*R.xopt + D;
        R.xopt = xopt;
        % solve overdetermined system to get multipliers for equalities
        % H*x + f + Am'*lambda_ineq + Ae'*lambda_eq = 0
        lambda_eq = zeros(S.me,1);
        lambda_eq(kept_rows.eq) = -Ae'\(S.H*R.xopt + S.f + Am'*lambda_ineq);
        
        % extend multipliers
        R.lambda.eqlin = lambda_eq;
    else
        R.lambda.eqlin = []; 
    end
            
    
else
    %% solve LCP directly
    if S.test
        [z,w,basis,exfl] = lcp(S.M, S.q);
    else        
        [z,w,basis,exfl] = lcp(S.M, S.q, MPTOPTIONS.modules.solvers.lcp );
    end

    R.xopt = z;
    R.lambda = w;
    
end

% this call of LCP solver is disabled, because the size of the problem is
% 2*n+m which is larger than with constraints factorization

%     if strcmpi(S.problem_type,'QP')
%         % if the problem is QP with H>0, then we can call LCP directly
%         v = eig(S.H);
%         if ~( any(abs(imag(v)) > MPTOPTIONS.abs_tol) || any(real(v) < -MPTOPTIONS.abs_tol) )
%             % Hessian is positive semidefinite
%             
%             % merge constraints
%             A = [-S.A; eye(S.n); -eye(S.n)];
%             b = [-S.b; S.lb; -S.ub];
%             
%             
%             if S.me==0
%                 % no equalities
%                 M = A*(S.H\A'); 
%                 q = -A*(S.H\S.f)-b;
%                 
%                 [z,w,basis,exfl] = lcp(full(M), full(q), MPTOPTIONS.modules.solvers.lcp );
%                 R.xopt = S.H\(A'*z - S.f);
%                 R.lamda = z(1:S.m);
%                 
%             else
%                 % equalities
%                 M = [A zeros(S.m+2*S.n,S.me)] * ( [S.H S.Ae'; -S.Ae zeros(S.me)] \ [A'; zeros(S.me,S.m+2*S.n)] );
%                 q = -[A zeros(S.m+2*S.n,S.me)] * ( [S.H S.Ae'; -S.Ae zeros(S.me)] \ [S.f; S.be] ) - b;
%                 
%                 [z,w,basis,exfl] = lcp(full(M), full(q), MPTOPTIONS.modules.solvers.lcp );
%                 xx = [S.H S.Ae'; -S.Ae zeros(S.me)] \ [A'*z - S.f; -S.be];
%                 R.xopt = xx(1:S.n);
%                 % merging multipliers on ineqality and equality constraints
%                 R.lambda = [z(1:S.m); xx(S.n+1:end)];
%                 
%             end
%             

% recalculate the objective function for LP, QP
R.obj = [];
if any(strcmpi(S.problem_type,{'LP','QP'}))
    if ~isempty(R.xopt) && ~isempty(S.H)
        R.obj = 0.5*R.xopt'*S.H*R.xopt + S.f'*R.xopt;
    elseif ~isempty(R.xopt) && ~isempty(S.f)
        R.obj = S.f'*R.xopt;
    end
end

switch exfl
    case 1
        R.how = 'ok';
        if S.test
            R.exitflag = 1;
        else
            R.exitflag = MPTOPTIONS.OK;
        end
    case -1
        R.how = 'infeasible';
        if S.test
            R.exitflag = 2;
        else
            R.exitflag = MPTOPTIONS.INFEASIBLE;
        end
    case -2
        R.how = 'unbounded';
        if S.test
            R.exitflag = 3;
        else
            R.exitflag = MPTOPTIONS.UNBOUNDED;
        end
    case -3
        R.how = 'preterminated (due to time limit or maximum pivot limit)';
        if S.test
            R.exitflag = -1;
        else
            R.exitflag = MPTOPTIONS.ERROR;
        end
    case -4
        R.how = 'other (numerical) error';
        if S.test
            R.exitflag = -1;
        else
            R.exitflag = MPTOPTIONS.ERROR;
        end
end

end
