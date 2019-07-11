classdef Opt < handle & matlab.mixin.Copyable
    % Encapsulates data and solutions for LP/QP/pLP/pQP/LCP/pLCP problems
    %
    % opt = Opt(param, value,...)
    % where param/value pairs define LP/QP/pLP/pQP variables
    % J(th) = min 0.5*u'*H*u + (pF*th+f)'*u + th'*Y*th + C*th + c
    %         s.t.  A*u <= b  + pB*th
    %               Ae*u = be + pE*th
    %               lb  <= u <= ub
    % or the LCP/pLCP variables
    % w - M*z = q + Q*th, w,z  >= 0, w'*z = 0
    %
    % opt = Opt(Polyheron P, param, value, ...)
    %  param/value pairs define cost, P defines constraints
    %
    
    properties (GetAccess=public, SetAccess=private)
        % LP/QP/pLP/pQP variables
        % J(th) = min 0.5*u'*H*u + (pF*th+f)'*u + th'*Y*th + C*th + c
        %         s.t.  A*u <= b  + pB*th
        %               Ae*u = be + pE*th
        %               lb  <= u <= ub
        %               Ath*th <= bth
        
        A  = []; b  = []; pB = [];
        Ae = []; be = []; pE = [];
        H  = []; f  = []; pF = [];
        lb = []; ub = [];
        Y = []; C = []; c = [];
        Ath = []; bth = [];
        
        % LCP variables
        % w - M*z = q + Q*th, w,z  >= 0, w'*z = 0
        M  = []; q  = []; Q  = [];
    end
    
    % public properties
    properties (GetAccess = public, SetAccess = public)
        Data  % any user-defined data
    end
    
    
    properties (SetAccess = private)
        n  = 0; % Problem dimension
        m  = 0; % Number of inequalities
        me = 0; % Number of equality constraints
        d  = 0; % Number of parameters
        
        % Problem types given as strings:
        % "LCP" - linear complementarity problem
        % "LP" - linear problem
        % "QP" - quadratic problem
        % "MILP" - mixed integer linear problem
        % "MIQP" - mixed integer quadratic problem

        problem_type = '';
        vartype = ''; % type of variables C-continuous, I-integer, B-binary, S-semicontinuous, N-semiinteger        
        isParametric = false;
        
        recover = []; % Mapping from solved problem to original problem data
        varOrder = [];
        Internal = [];
	end
	
	properties
		solver = ''
	end
    
	methods
		function set.solver(obj, new_solver)
			% Opt.solver setter
			
			if ~ischar(new_solver)
				error('The solver must be a string.');
			end
			obj.solver = upper(new_solver);
		end
	end
	
    methods(Access = public)
        
        % Constructor
        function opt = Opt(varargin)
            if nargin > 0
                if isa(varargin{1},'lmi') || isa(varargin{1},'constraint')
                    % convert from YALMIP data
                    opt = opt.setYalmipData(varargin{:});
                elseif isa(varargin{1},'struct')
                    % convert from MPT2.6 data
                    if isfield(varargin{1},'G') && isfield(varargin{1},'W') && isfield(varargin{1},'E')
                        opt = opt.setMPT26Data(varargin{1});
                    else
                        % call normal constructor
                        opt = opt.setData(varargin{:});
                    end
                else
                    % call normal constructor
                    opt = opt.setData(varargin{:});
                end
                
                % validate
                opt.validate;
                    
            else
                error('Empty problems not supported.');
            end
        end
        
        function K = feasibleSet(obj, arg)
            % Computes the feasible set of a given parametric problem
            %
            % For a parametric problem
            %    min  J(z, x)
            %    s.t. A*z <= b + pB*x
            %         Ae*z = be + pE*x
            % the feasible set K is the polyhedron
            %   K = { x | \exists z s.t. A*z<=b+pB*x, Ae*z=be+pE*x }
            %
            % This method implements two procedures to compute K:
            %   1) if K=prob.feasibleSet() is called, the feasible set is
            %      calculated by projection (can be expensive)
            %   2) if K=prob.feasibleSet(regions) is called with "regions"
            %      being the critical regions of the parametric solution,
            %      then K is constructed as follows:
            %         For each facet of each region do:
            %          a) compute the center of the facet
            %          b) take a small step accross the facet
            %          c) solve the problem for the new point
            %          d) if the problem is infeasible, add the facet to
            %             the feasible set 
            %
            % Syntax:
            %   K = prob.feasibleSet()
            %   K = prob.feasibleSet(method)
            %   K = prob.feasibleSet(regions)
            %
            % Inputs:
            %      prob: parametric problem as an Opt object
            %   regions: (optional) critical regions of the parametric
            %            solution
            %    method: (optional) string identificator of the projection
            %            method to use (see help Polyhedron/projection). By
            %            default we use the 'ifourier' method.
            %
            % Output:
            %         K: feasible set as a redundant H-polyhedron
            
            global MPTOPTIONS
            
            narginchk(1, 2);
            if ~obj.isParametric
                error('The problem is not parametric.');
            end
            if nargin==2
                if isa(arg, 'Polyhedron')
                    use_projection = false;
                    regions = arg;
                elseif isa(arg, 'char')
                    use_projection = true;
                    method = arg;
                else
                    error('The input must be either a string or an array of Polyhedron objects.');
                end
            else
                % default projection method
                use_projection = true;
                method = 'ifourier';
            end
            
            if use_projection
                % compute the feasible set via projection
                
                if isequal(lower(obj.problem_type), 'lcp') && ...
                        ~isfield(obj.Internal, 'constraints')
                    % This is a pLCP that was defined manually
                    %
                    % LCP constraints:
                    %   I*w - M*z = q + Q*x
                    %   w >= 0, z >= 0
                    %   Ath*x <= bth
                    
                    [nw, nz] = size(obj.M);
                    nb = length(obj.bth);
                    % y = [x; w; z]
                    Ae = [-obj.Q, eye(nw), -obj.M];
                    be = obj.q;
                    A = [obj.Ath, zeros(nb, nw+nz); ...
                        zeros(nw, obj.d), -eye(nw), zeros(nw, nz); ...
                        zeros(nz, obj.d+nw), -eye(nz)];
                    b = [obj.bth; zeros(nw+nz, 1)];
                    ZX = Polyhedron('Ae', Ae, 'be', be, 'A', A, 'b', b);

                else
                    % LP/QP constraints:
                    %     A*z <= b + pB*x
                    %    Ae*z == be + pE*x
                    %   Ath*x <= bth
                    
                    if isfield(obj.Internal, 'constraints')
                        % This is a pLCP that was generated by
                        % Opt/qp2lcp(). Convert it to the original pQP
                        % formulation to reduce the number of variables
                        obj = Opt(obj.Internal.constraints);
                    end
                    
                    if obj.me>0
                        % eliminate equalities, but do it on a copy of the object
                        obj = obj.copy();
                        obj.eliminateEquations();
                    end
                    
                    % construct the polyhedron 
                    % { (x, z) | A*x<=b+pB*x, Ath*x<=bth }
                    ZX = Polyhedron([-obj.pB, obj.A; ...
                        obj.Ath, zeros(length(obj.bth), obj.n)], ...
                        [obj.b; obj.bth]);
                end
                
                if ZX.isEmptySet()
                    % the feasible set is empty
                    K = Polyhedron.emptySet(obj.d);
                    return
                end
                
                % the feasible set is given as the projection of ZX onto
                % the parametric space
                K = ZX.projection(1:obj.d, method);
                
            else
                % construct the feasible set via critical regions
                
                % length of step over the facet
                step_size = MPTOPTIONS.rel_tol*10;
                Hf = [];
                t = tic;
                n_fails = 0;
                for i = 1:length(regions)
                    % for each region
                    if toc(t) > MPTOPTIONS.report_period
                        fprintf('progress: %d/%d\n', i, length(regions));
                        t = tic;
                    end
                    % make sure we have the minimal H-representation
                    % (we do redundancy elimination here since it can be
                    % costly in high dimensions; hence we want to give the
                    % uer an appropriate progress report)
                    regions(i).minHRep();
                    for j = 1:length(regions(i).b)
                        % for each facet of the i-th region:
                        % 1) compute a point on the facet
                        lpsol = regions(i).chebyCenter(j);
                        if lpsol.exitflag == MPTOPTIONS.OK
                            % 2) compute the point accross the j-th facet
                            x = lpsol.x + regions(i).A(j, :)'/norm(regions(i).A(j, :)')*step_size;
                            % 3) and solve the problem for the new point
                            qpsol = obj.solve(x);
                            if qpsol.exitflag ~= MPTOPTIONS.OK
                                % 4) infeasible => add this facet to the feasible set
                                Hf = [Hf; regions(i).H(j, :)];
                            end
                        else
                            % numerical problems
                            n_fails = n_fails + 1;
                        end
                    end
                end
                if n_fails > 0
                    fprintf('WARNING: failed to compute points on %d facet(s)\n', n_fails);
                end
                if isempty(Hf)
                    % numerical problems, return R^n
                    K = Polyhedron.fullSpace(regions(1).Dim);
                else
                    K = Polyhedron(Hf(:, 1:end-1), Hf(:, end));
                end
            end
        end
        
        function obj = minHRep(obj)
            % Removes redundant inequalities
            %
            % Given a parametric optimization problem with constraints
            %      A*z <= b + pB*theta
            % this method removes the redundant inequalities.
            
            assert(obj.d>0, 'The problem is not parametric.');
            
            P = Polyhedron([obj.A, -obj.pB], obj.b);
            P.minHRep();
            obj.A = P.A(:, 1:obj.n);
            obj.pB = -P.A(:, obj.n+1:end);
            obj.b = P.b;
            obj.m = length(obj.b);
        end
        
        function CR = getRegion(obj, A, implicit_string)
            % Constructs a critical region, the optimizers, and the cost function
            %
            % Syntax:
            %   CR = pqp.getRegion(A)
            %   CR = pqp.getRegion(A, 'implicit') - for implicit representation
            %
            % Inputs:
            %   pqp: pqp problem as an instance of the Opt class
            %     A: incides of the active inequalities
            %
            % Output:
            %     R: critical region as a Polyhedron object with functions "primal"
            %        and "obj"

            narginchk(2, 3);
            assert(isequal(obj.problem_type, 'QP'), 'Only pQPs are supported.');
            assert(obj.d>0, 'The problem must have parameters.');
            assert(obj.me==0, 'Equalities are not supported.');
            if nargin<3
                implicit_string = '';
            end
            implicit = ~isempty(implicit_string);
            
            A(A==0) = []; % depad zeros
            
            % TODO: resolve primal degeneracy
            assert(numel(A)<=obj.n, 'Primally degenerate active set.');
            
            % special notion for no active constraints
            if isequal(A, 0)
                A = [];
            end
            
            % extract active sets
            Ga = obj.A(A, :);
            Ea = obj.pB(A, :);
            wa = obj.b(A, :);
            
            if rank(Ga)~=numel(A)
                error('Primally degenerate active set.');
            end
            
            % extract non-active sets
            all = 1:obj.m;
            N  = all(~ismembc(all, A));
            Gn = obj.A(N, :);
            En = obj.pB(N, :);
            wn = obj.b(N, :);
            
            % optimal dual variables
            alpha_1 = -Ga/obj.H*obj.pF - Ea;
            alpha_2 = -Ga/obj.H*Ga';
            beta    = -Ga/obj.H*obj.f - wa;
            
            % dual variables are an affine function of the parameters
            alpha_L = -alpha_2\alpha_1;
            beta_L  = -alpha_2\beta;
            
            % optimizer is an affine function of the parameters
            alpha_x = -obj.H\obj.pF - obj.H\Ga'*alpha_L;
            beta_x  = -obj.H\obj.f - obj.H\Ga'*beta_L;

            % cost function
            % TODO: check this
            Jquad = 0.5*alpha_x'*obj.H*alpha_x + obj.pF'*alpha_x + obj.Y;
            Jaff = beta_x'*obj.H*alpha_x + beta_x'*obj.pF + obj.f'*alpha_x + obj.C;
            Jconst = 0.5*beta_x'*obj.H*beta_x + obj.f'*beta_x + obj.c;
            J = QuadFunction(Jquad, Jaff, Jconst);

            % Critical region
            crH = [Gn*alpha_x - En; -alpha_L];
            crh = [-Gn*beta_x + wn; beta_L];
            
            if ~implicit
                % include Ath*x<=bth bounds
                crH = [crH; obj.Ath];
                crh = [crh; obj.bth];
                CR = Polyhedron(crH, crh);
                CR.Data.ActiveSet = A;
            else
                data.ActiveSet = A;
                data.OptProb = obj;
                data.Primal.F = alpha_x;
                data.Primal.g = beta_x;
                data.DualIneq.F = alpha_L;
                data.DualIneq.g = beta_L;
                CR = IPDPolyhedron(data);
            end

            % Project primal optimizer back on equalities
            if ~isempty(obj.recover)
                Lprimal = obj.recover.Y*[alpha_x, beta_x] + obj.recover.th;
                alpha_x = Lprimal(:, 1:end-1);
                beta_x = Lprimal(:, end);
            end
            % Extract only requested variables from the primal optimizer
            if ~isempty(obj.varOrder) && ~isempty(obj.varOrder.requested_variables)
                alpha_x = alpha_x(obj.varOrder.requested_variables, :);
                beta_x = beta_x(obj.varOrder.requested_variables);
            end
            z = AffFunction(alpha_x, beta_x);
            % lagrange multipliers: aL*x+bL
            aL = zeros(size(obj.pB));
            bL = zeros(size(obj.b));
            aL(A, :) = alpha_L;
            bL(A) = beta_L;
            L = AffFunction(aL, bL);
            
            % attach functions
            CR.addFunction(z, 'primal');
            CR.addFunction(J, 'obj');
            CR.addFunction(L, 'dual-ineqlin');
        end
        
        function [A, sol] = getActiveSetForPoint(obj, theta)
            % Returns constraints that are active for a given parameter
            %
            % [A, SOL] = prob.getActiveSetForPoint(THETA) returns the
            % indices A of inequality constraints that are active at the
            % optimum for the given parameter THETA. The SOL output
            % determines feasibility in SOL.how and the primal optimizer
            % in SOL.xopt.
            
            global MPTOPTIONS
            assert(~isequal(obj.problem_type, 'LCP'), 'LCP problems are not supported.');
            assert(obj.isParametric, 'The problem must have parameters.');
            assert(isequal(size(theta), [obj.d 1]), 'The parameter must be %dx1.', obj.d);
            sol = obj.solve(theta);
            if sol.exitflag==MPTOPTIONS.OK
                A = find(abs(obj.A*sol.xopt-obj.b-obj.pB*theta)<MPTOPTIONS.abs_tol);
            else
                % infeasible
                A = NaN;
            end
        end
        
        function [result, nlps] = checkActiveSet(obj, A)
            % Checks optimality/fesibility of a given active set W
            %
            % pQP formulation:
            %
            %    min_z 0.5*z'*H*z + (pF*x+f)'*z + x'*Y*x + C'*x + c
            %     s.t.   A*z <=  b + pB*x
            %           Ae*z == be + pE*x
            %          Ath*x <= bth
            %             lb <= z <= ub
            %
            % KKT conditions:
            %         stanionarity: H*z + pF*x + f + Aa'*La = 0
            %     compl. slackness: Aa*z - ba - pBa*x = 0
            %   primal feasibility: An*z - bn - pBn*x < 0
            %                       Ae*z - be - pE*x  = 0
            %     dual feasibility: La >= 0
            %
            % where Aa, ba, pBa are matrices formed by taking rows indexed by W from
            % the corresponding matrix; An, bn, pBn contain only the inactive rows.
            %
            % To certify that W is an optimal active set, we solve an LP:
            %
            %    max  t
            %    s.t. H*z + pF*x + f + Aa'*La  =  0  (optimality)
            %                     Aa*z - pBa*x = ba (complementarity slackness)
            %                     Ae*z -  pE*x = be (primal feasibility, equalities)
            %                An*z - pBn*x + t <= bn  (primal feasibility, inequalities)
            %                              La >= t   (dual feasibility)
            %
            % * if the LP is feasible with t>0, the active set is optimal,
            %   and the dual variables satisfy the strict complementarity
            %   slackness condition (i.e., La>0) 
            % * if the LP is feasible with t=0, the active set if optimal,
            %   but the dual variables violate the strict complementarity
            %   slackness condition (i.e., La>=0)
            % * if the LP is feasible without the optimality constraints,
            %   the active set is feasible but not optimal 
            %   (not checked if card(W)=nz)
            % * otherwise the active set is infeasible
            % * t>0 is checked as t>=-MPTOPTIONS.zero_tol
            %
            % Before solving the LP we first check whether A(W, :) has full
            % row rank. If not, the active set violates the linear
            % independence constraint qualification condition, i.e., some
            % rows of A(W, :) are linearly dependent. In that case we exit
            % immediately with result=-2
            %
            % Usage:
            %   [status, nlps] = obj.checkActiveSet(W)
            %
            % Inputs:
            %   obj: matrices of the (p)QP problem as an Opt object
            %     W: active set to investigate (=indices of active constraints)
            %
            % Outputs:
            %   result: flag indicating the status of the active set
            %          2: optimal, strict complementary slackness holds
            %          1: optimal, strict complementary slackness violated
            %          0: feasible, but not optimal
            %         -1: infeasible
            %         -2: rank defficient (LICQ violation)
            %         -3: undecided
            %       nlps: number of LPs solved
            
            % TODO:
            %   3: optimal, primal full dim, L>0
            %   2: optimal, primal full dim, L>=0
            %   1: optimal, primal low dim, L>=0
            %   0: feasible
            %  -1: infeasible
            %  -2: rank defficient
            %  -3: undecided
            
            global MPTOPTIONS
            
            % TODO: support pLPs, pLCPs, QPs, LPs, LCPs
            assert(isequal(lower(obj.problem_type), 'qp'), 'Only (p)QPs are supported for now.');
            
            A(A==0) = []; % depad zeros
            
            nlps = 0;
            Ga = obj.A(A, :);
            % check rank of Ga first
            rGa = rank(Ga'*Ga, MPTOPTIONS.abs_tol);
            if rGa < length(A)
                % Ga is rank deficient
                result = -2;
                return
            end
            
            % determine the index set of inactive constraints
            all = 1:obj.m;
            %N = setdiff(1:obj.ni, A);
            N = all(~ismembc(all, A)); % much faster version of setdiff for sorted arrays
            
            Gn = obj.A(N, :);
            Sa = obj.pB(A, :);
            Sn = obj.pB(N, :);
            wa = obj.b(A, :);
            wn = obj.b(N, :);
            nA = length(A);
            nN = length(N);
            
            % optimization variables: [t; z; x; La]
            
            % max t
            lp.f = [-1 zeros(1, obj.n+obj.d+nA)];
            
            % equality constraints:
            %   H*z + pF*x + f + Ga'*La =  0  (optimality)
            %               Ga*z - Sa*x = wa  (complementarity slackness)
            %               Ge*z - Se*x = we  (primal feasibility, equalities)
            lp.Ae = [ zeros(obj.n, 1), obj.H, obj.pF, Ga'; ...
                zeros(nA, 1), Ga, -Sa, zeros(nA); ...
                zeros(obj.me, 1), obj.Ae, -obj.pE, zeros(obj.me, nA)];
            lp.be = [ -obj.f; wa; obj.be ];
            
            % inequality constraints:
            %   Gn*z - Sn*x + t  <= wn  (primal feasibility, inequalities)
            %                 La >= t   (dual feasibility)
            %              Ath*x <= bth
            nb = length(obj.bth);
            lp.A = [ ones(nN, 1), Gn, -Sn, zeros(nN, nA); ...
                ones(nA, 1), zeros(nA, obj.n), zeros(nA, obj.d), -eye(nA); ...
                zeros(nb, 1+obj.n), obj.Ath, zeros(nb, nA)];
            
            % TODO: strictly speaking, we should not enforce La>=t
            %lp.A = [ ones(nN, 1), Gn, -Sn, zeros(nN, nA); ...
            %    zeros(nA, 1), zeros(nA, obj.n), zeros(nA, obj.d), -eye(nA); ...
            %    zeros(nb, 1+obj.n), obj.Ath, zeros(nb, nA)];
            
            lp.b = [ wn; zeros(nA, 1); obj.bth];
            
            % lower/upper bounds on [t; z; x; La]
            % (includes lb <= z <= ub)
            lp.lb = [-Inf; obj.lb; -Inf(obj.d, 1); zeros(nA, 1)];
            lp.ub = [Inf; obj.ub; Inf(obj.d, 1); Inf(nA, 1)];
            
            % avoid sanity checks in mpt_solve()
            lp.quicklp = true;
            
            nlps = nlps + 1;
            sol = mpt_solve(lp);
            if isequal(sol.how, 'ok')
                if sol.xopt(1) > MPTOPTIONS.rel_tol
                    % feasible, optimal, LICQ holds
                    result = 2;
                elseif sol.xopt(1) > -MPTOPTIONS.zero_tol
                    % feasible, optimal, LICQ does not hold
                    result = 1;
                else
                    % undecided
                    result = -3;
                end
            else
                % undecided
                result = -3;
            end
            
            if (result==-3 || result==-1) && numel(A)<obj.n
                % solve the LP without optimality conditions (not necessary if we are
                % on the last level)
                
                % TODO: solve just one LP with soft equality optimality constraint
                %       H*z + pF*x + f + Ga'*La = s, min Qs*norm(s, 1)
                lp.Ae = lp.Ae(obj.n+1:end, :);
                lp.be = lp.be(obj.n+1:end);
                nlps = nlps + 1;
                sol = mpt_solve(lp);
                if isequal(sol.how, 'ok') && sol.xopt(1)>-MPTOPTIONS.zero_tol
                    % feasible, but not optimal
                    result = 0;
                    
                    if false
                        % EXPERIMENTAL CODE
                        %
                        % how far away are we from optimality?
                        q = 0;
                        z=sdpvar(obj.nz, 1);
                        x=sdpvar(obj.nx, 1);
                        La=sdpvar(nA, 1);
                        t=sdpvar(1, 1);
                        e=sdpvar(1, 1);
                        
                        % min e - q*t
                        %
                        % -e <= H*z + Ga'*La <= e (weak optimality)
                        % Ga*z - Sa*x - wa =   0  (complementarity slackness)
                        % Gn*z - Sn*x + t <= wn   (primal feasibility)
                        % La >= t                 (dual feasibility)
                        % t <= 1
                        F = [ -e <= obj.H*z+Ga'*La <= e; ...
                            Ga*z - Sa*x - wa == 0; ...
                            Gn*z - Sn*x + t <= wn; ...
                            La >= t; ...
                            0 <= t <= 1; e >= 0 ];
                        obj = e-q*t;
                        info = solvesdp(F, obj, sdpsettings('verbose', 0));
                        z=double(z);
                        x=double(x);
                        La=double(La);
                        t=double(t);
                        e=double(e);
                        if e>500
                            result = -1;
                        end
                        %OPTGAP{end+1} = struct('A', A, 'e', e, 't', t);
                    end
                else
                    % infeasible
                    result = -1;
                end
            end
        end
        
        function [Aopt, Adeg, Afeas, Ainfeas, nlps] = enumerateActiveSets(obj, varargin)
            % Enumerates all optimal combinations of active sets
            %
            %   [Aopt, Adeg, Afeas, Ainfeas, nLPs] = obj.enumerateActiveSets()
            %
            % Aopt: optimal active sets
            % Adeg: degenerate active sets
            % Afeas: feasible active sets
            % Ainfeas: infeasible active sets
            % nLPs: number of LPs needed to find the sets
            %
            %    obj.enumerateActiveSets('opt1', value1, 'opt2', value2, ...)
            %
            % Supported options:
            %    'verbose': level of verbosity
            %    'prune_infeasible': if true (default), infeasible active
            %                        sets are not followed
            %    'exclude': list of indices of constraints to exclude
            
            global MPTOPTIONS
            if isempty(MPTOPTIONS)
                MPTOPTIONS=mptopt;
            end
            p = inputParser;
            p.addParamValue('verbose', MPTOPTIONS.verbose);
            p.addParamValue('report_period', MPTOPTIONS.report_period);
            p.addParamValue('prune_infeasible', true);
            p.addParamValue('exclude', []);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            options = p.Results;

            % TODO: support pLPs, pLCPs, LPs, LCPs
            assert(isequal(lower(obj.problem_type), 'qp'), 'Only (p)QPs are supported for now.');
            
            start_t = clock;
            
            % is the case with no active constraints optimal?
            [feasible, nlps] = obj.checkActiveSet([]);
            if feasible>=2
                % yep, include the case into list of non-degenerate optimal active sets
                Aopt = {0};
                Afeas = {0};
            else
                % nope, start with an empty list
                Aopt = {[]};
                Afeas = {-1};
            end
            Adeg = {[]};
            Ainfeas = {[]};
            
            if options.verbose>=0
                fprintf('Level    Total Candidates Optimal Degenerate Feasible Infeasible Rankdef       LPs\n');
            end
            
            % since we have "nz" optimization variables and a strictly convex QP, at
            % most "nz" constraints will be active at the optimum
            for i = 1:obj.n
                if options.verbose>=0
                    fprintf('%2d/%2d', i, obj.n);
                end
                % explore all nodes in this level, provide unique list of feasible
                % constraints and all previously discovered infeasible combinations
                [Ao, Ad, Af, Ai, nlp] = sub_exploreLevel(obj, i, Afeas{end}, Ainfeas, options);
                nlps = nlps + nlp;
                Aopt{i+1} = Ao;    % optimal active sets
                Adeg{i+1} = Ad;    % degenerate active sets
                Afeas{i+1} = Af;   % feasible active sets
                Ainfeas{i+1} = Ai; % infeasible active sets
            end
            if options.verbose>=0
                fprintf('...done (%.1f seconds, %d LPs)\n', etime(clock, start_t), nlps);
            end
            
            % convert cells to zero-padded matrices
            Aopt = sub_cell2mat_pad(Aopt, obj.n);
            Adeg = sub_cell2mat_pad(Adeg, obj.n);
            Afeas = sub_cell2mat_pad(Afeas, obj.n);
            Ainfeas = sub_cell2mat_pad(Ainfeas, obj.n);
            if Afeas(1, 1)==-1
                % remove a dummy row
                Afeas = Afeas(2:end, :);
            end
        end

    end
    
%     methods
%         %% SET methods
%         % check if vartype is correct
%         function set.vartype(obj,val)
%             if ~isempty(val)
%                 if isnumeric(val)
%                     % convert to char if it is numeric
%                     val = char(val);
%                 end
%                 if ~isvector(val) || ~ischar(val)
%                     error('The argument must be a vector of strings.');
%                 end
%                 % checking if string is correct
%                 for i=1:length(val)
%                     if ~any(strcmpi(val(i),{'C','I','B','S','N'}))
%                         %C-continuous, I-integer, B-binary, S-semicontinuous, N-semiinteger
%                         error('Given string does not belong to gropu "C-continuous, I-integer, B-binary, S-semicontinuouos, N-semiinteger.');
%                     end
%                 end
%                 
%                 obj.vartype = val;
%             end            
%         end
%         
%     end
    
    methods (Access=protected)
        function U = copyElement(obj)
            % Creates a copy of the union
            %
            %   copy = U.copy()
            
            % Note: matlab.mixin.Copyable.copy() automatically properly
            % copies arrays and empty arrays, no need to do it here.
            % Moreover, it automatically creates the copy of proper class.
            U = copyElement@matlab.mixin.Copyable(obj);

            % we don't know what type of arguments can be put in the future
            % to Internal properties, so we check if a field contains a
            % "copy" method
            if isstruct(obj.Internal)
                nf = fieldnames(obj.Internal);
                for i=1:numel(nf)
                    x = obj.Internal.(nf{i});
                    if isobject(x) && ismethod(x, 'copy');
                        U.Internal.(nf{i}) = x.copy();
                    end
                end
            else
                if isobject(obj.Internal) && ismethod(obj.Internal,'copy');
                    U.Internal = obj.Internal.copy();
                end
            end
            % we don't know what type of arguments can be stored inside
            % Data, so we check if it contains a "copy" method - then use
            % it to create new object.
            if isstruct(obj.Data)
                nd = fieldnames(obj.Data);
                for i=1:numel(nd)
                    x = obj.Data.(nd{i});
                    if isobject(x) && ismethod(x,'copy');
                        U.Data.(nd{i}) = x.copy();
                    end
                end
            else
                if isobject(obj.Data) && ismethod(obj.Data,'copy');
                    U.Data = obj.Data.copy;
                end                
            end
            
        end
 
    end
end

function M = sub_cell2mat_pad(A, n)

n_elements = 0;
for i = 1:numel(A)
    n_elements = n_elements + size(A{i}, 1);
end
M = zeros(n_elements, n);
start_idx = 1;
for i = 1:numel(A)
    end_idx = start_idx+size(A{i}, 1);
    M(start_idx:end_idx-1, 1:size(A{i}, 2)) = A{i};
    start_idx = end_idx;
end

end

function [Aopt, Adeg, Afeasible, Ainfeasible, nlps] = sub_exploreLevel(pqp, level, feasible, infeasible, options)
% Checks feasibility/optimality of each n-combination of active constraints
%
% Syntax:
% -------
%
% [Ao, Ad, Af, Ai] = exploreLevel(pqp, level, feasible, infeasible, options)
%
% Inputs:
% -------
%             pqp: matrices of the pQP formulation
%           level: index of the level (integer)
%        feasible: cell array of feasible active sets at the previous level
%      infeasible: cell array of infeasible active sets at each level
% options.verbose: if >=0, progress will be displayed
%
% Outputs:
% --------
%           Ao: m-by-n matrix of optimal active sets without SCS violation (n=level)
%           Ad: m-by-n matrix of optimal active sets with SCS violation
%           Af: m-by-n matrix of non-optimal feasible active sets
%           Ai: m-by-n matrix of infeasible active sets
%               (includes rank-defficient active sets)
%         nlps: number of LPs solved on this level

AllFeasible = unique(feasible(:));
% pre-alocate arrays
Aopt = zeros(size(feasible, 1)*length(AllFeasible), level);
Adeg = zeros(size(feasible, 1)*length(AllFeasible), level);
Afeasible = zeros(size(feasible, 1)*length(AllFeasible), level);
Ainfeasible = zeros(size(feasible, 1)*length(AllFeasible), level);
n_feasible = 0;
n_opt = 0;
n_deg = 0;
n_infeasible = 0;
n_rankdef = 0;
n_candidates = 0;
nlps = 0;
n_pruned = 0;
t=tic;
first_display = true;
for i = 1:size(feasible, 1)
    if level==1
        candidates = 1:pqp.m;
    else
        candidates = AllFeasible(AllFeasible>max(feasible(i, :)));
    end
    if ~isempty(options.exclude)
        % remove user-defined constraints to be excluded
        candidates = setdiff(candidates, options.exclude);
    end
    for j = candidates(:)'
        if level==1
            Atry = j;
        else
            Atry = [feasible(i, :), j];
        end
        
        % remove previously known infeasible combinations
        %
        % we only need to do this from the 3rd level upwards, since on the 2nd
        % level the lists of feasible and infeasible constraints are mutually
        % exclusive, hence nodes are not poluted by any infeasible constraints
        do_check = true;
        if level>2 && options.prune_infeasible
            % check if Atry contains any sequence listed in "infeasible"
            for j = 2:length(infeasible)
                m = ismembc(infeasible{j}, Atry); % fastest implementation
                if any(sum(m, 2)==size(infeasible{j}, 2))
                    % Atry contains all entries from at least one row of
                    % "infeasible", therefore it can be removed
                    do_check = false;
                    n_pruned = n_pruned+1;
                    break
                end
            end
        end
        if do_check
            n_candidates = n_candidates + 1;
            [result, nlp] = pqp.checkActiveSet(Atry);
            nlps = nlps+nlp;
            %  2: optimal, no SCS violation
            %  1: optimal, SCS violated
            %  0: feasible, but not optimal
            % -1: infeasible
            % -2: rank defficient, LICQ violated
            % -3: undecided
            % split active sets into feasible/infeasible/optimal/degenerate
            if result>=0
                n_feasible = n_feasible+1;
                Afeasible(n_feasible, :) = Atry;
                if result==1
                    n_deg = n_deg + 1;
                    Adeg(n_deg, :) = Atry;
                elseif result>=2
                    n_opt = n_opt + 1;
                    Aopt(n_opt, :) = Atry;
                end
            end
            if result<0
                n_infeasible = n_infeasible + 1;
                Ainfeasible(n_infeasible, :) = Atry;
                if result==-2
                    n_rankdef = n_rankdef + 1;
                end
            end
        end
    end
    if options.verbose>=0 && toc(t)>options.report_period
        if ~first_display
            fprintf(repmat('\b', 1, 9));
        end
        fprintf('%8d%%', min(100, ceil(100*i/size(feasible, 1))));
        t=tic;
        first_display = false;
    end
end
Afeasible = Afeasible(1:n_feasible, :);
Ainfeasible = Ainfeasible(1:n_infeasible, :);
Aopt = Aopt(1:n_opt, :);
Adeg = Adeg(1:n_deg, :);

if options.verbose>=0
    % delete progress meter
    if ~first_display
        fprintf(repmat('\b', 1, 9));
    end
    % display progress
    fprintf('%9d   %8d', n_candidates+n_pruned, n_candidates);
    fprintf('%8d   %8d %8d   %8d%8d  %8d\n', ...
        n_opt, n_deg, n_feasible-n_opt-n_deg, n_infeasible-n_rankdef, ...
        n_rankdef, nlps);
end

end
