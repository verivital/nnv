function ret = mpt_call_mpqp(S)
%
% a gateway routine to parametric QP solver of MPT2.6
%

global MPTOPTIONS
if isempty(MPTOPTIONS)
    MPTOPTIONS = mptopt;
end

if ~isa(S,'Opt')
    error('mpt_call_mpqp: Input argument must be an "Opt" class.');
end

if ~strcmpi(S.problem_type,'QP')
    error('mpt_call_mpqp: MPQP solver does not solve %s problems!',S.problem_type);
end


% store the original problem
opt = S.copy;
% make a copy of the original problem before we modify it via
% eliminateEquations()
S = S.copy();

if S.me>0
    % eliminate equations first
    S.eliminateEquations;
end

% In validation of Opt class there are prepreprocessing functions
% for tightening the bounds on the parametric solution. One of the task performs
% extraction of lower and upper bounds on the variables from equation 
% G*U <= W + E*th and puts them into separate fields. We need to put these 
% bounds back.

Matrices.G = S.A;
Matrices.W = S.b;
Matrices.E = S.pB;
ilb = (S.lb==-Inf) | (S.lb<-MPTOPTIONS.infbound);
iub = (S.ub==Inf)  | (S.ub>MPTOPTIONS.infbound);
if any(~ilb)
    % put ones at the positions where there is ub
    L = -eye(S.n);
    Matrices.G = [Matrices.G; L(~ilb,:)];
    Matrices.W = [Matrices.W; -S.lb(~ilb)];
    Matrices.E = [Matrices.E; zeros(nnz(~ilb),S.d)];
end
if any(~iub)
    % put ones at the positions where there is ub
    U = eye(S.n);
    Matrices.G = [Matrices.G; U(~iub,:)];
    Matrices.W = [Matrices.W; S.ub(~iub)];
    Matrices.E = [Matrices.E; zeros(nnz(~iub),S.d)];
end
Matrices.H = S.H;
Matrices.F = S.pF';
Matrices.Cf = S.f';
Matrices.Y = S.Y;
Matrices.Cx = S.C;
Matrices.Cc = S.c;
Matrices.bndA = S.Ath;
Matrices.bndb = S.bth;

if MPTOPTIONS.verbose >= 1
    disp('Calling mpt_mpqp_26 with default options...')
end
start_time = clock;
[r.Pn,r.Fi,r.Gi,r.activeConstraints,r.Phard,r.details]=mpt_mpqp_26(Matrices);

% Re-order variables if this came from YALMIP
P = speye(opt.n);
if ~isempty(opt.varOrder)
    P = P(opt.varOrder.requested_variables,:);
end

% convert @polytope to @Polyhedron
reg = toPolyhedron(r.Pn);

for i=1:length(reg)
    
    % add only not empty regions (with region_tol)
    if ~reg(i).isEmptySet
        %regions.add(reg);
        
        % y = Fi*th + Gi = [Fi Gi]*[th;1]
        % x = Y*y + T*[th;1];
        % x = (Y*[Fi Gi]+T)*[th;1]
        
        if opt.me
            % Compute affine mapping from parameter to primal        
            T = S.recover.Y*[r.Fi{i} r.Gi{i}] + S.recover.th;
        else
            T = [r.Fi{i} r.Gi{i}];
        end
        
        % compute primal variables
        Lprimal = P*T;
        reg(i).addFunction(AffFunction(Lprimal(:,1:end-1),Lprimal(:,end)),'primal');
        
        % compute the objective value
        Y = T(:,1:end-1);
        R = T(:,end);

        if ~isempty(opt.H)
            qt = 0.5*Y'*opt.H*Y + opt.pF'*Y + opt.Y;
            lt = R'*opt.H*Y + R'*opt.pF + opt.f'*Y + opt.C;
            at = 0.5*R'*opt.H*R + opt.f'*R + opt.c;
            reg(i).addFunction(QuadFunction(qt,lt,at),'obj');
        else
            lt = R'*opt.pF + opt.f'*Y + opt.C;
            at = opt.f'*R + opt.c;
            reg(i).addFunction(AffFunction(lt,at),'obj');
        end

    end
end

ret.xopt = PolyUnion('Set',reg,'Domain', toPolyhedron(r.Phard),...
	'Convex',true,'Overlaps',false,'Bounded',true,'Fulldim',true,'Connected',true);
ret.xopt.setInternal('convexHull', toPolyhedron(r.Phard));
ret.mpqpsol = r;
if ret.xopt.Num>0
	ret.exitflag = MPTOPTIONS.OK;
	ret.how = 'ok';
else
	ret.exitflag = MPTOPTIONS.INFEASIBLE;
	ret.how = 'infeasible';
end
ret.stats.solveTime = etime(clock, start_time);

end
