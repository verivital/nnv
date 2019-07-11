function options = mpt_solvers_options
%
% this function returns initial settings for solvers in this directory
%
% User is allowed to change these settings later using "mptopt" class.

%% CLP options
options.clp.solver = 1; %  [1 (primal)| 2 (dual), (default 1)].
options.clp.maxnumiterations = 99999999; % [int>=0 (default 99999999)]
options.clp.maxnumseconds = 3600; % [int>=0 (default 3600)]
options.clp.primaltolerance  = 1e-7; % [double>=0 (default 1e-7)]
options.clp.dualtolerance    = 1e-7; % [double>=0 (default 1e-7)]
options.clp.primalpivot = 1; %  [1 (steepest) | 2 (Dantzig) (default 1)]
options.clp.dualpivot = 1; % [1 (steepest) | 2 (Dantzig) (default 1)]
options.clp.verbose = 0; % [0|1|... (default 0)]

%% CPLEX options
% cplexint interface
options.cplexint.verbose = 0; % [0|1|2 (default 0)]
options.cplexint.logfile = 0; % [0|1 (default 0)]
options.cplexint.lic_rel = 1e3; % after how many runs to release the license

% cplex interface by IBM 
if exist('Cplex','file')==6
    matlab_version = version('-release');
    if str2double(matlab_version(1:4))>=2016
        % for CPLEX 12.7
        options.cplex = cplexoptimset;
    else
        options.cplex = cplexoptimset('cplex');
    end
    options.cplex.Display = 'off';
    c=Cplex;
    if strcmp(c.Version,'12.2.0.0')
        % for version 12.2 we need to add this option to disable verbosity,
        % but the selection of LP/QP method is not accepted which gives
        % wrong results in some tests
        options.cplex.Diagnostics = 'off';
        disp('Please, update your CPLEX software for version 12.4 or higher because this version does not accept modified options.');
    end
end
% parsing the all set of options slows down the computation
% tremendously!!!, but if we do not provide values to all options, then the
% routine does not accept the new values

% 0 Automatic
% 1 Primal Simplex
% 2 Dual Simplex
% 3 Network Simplex
% 4 Barrier
% 5 Sifting
% 6 Concurrent
options.cplex.lpmethod = 2; % dual-simplex
options.cplex.qpmethod = 2; % dual-simplex

%% PLCP options
options.plcp.bfs = true;
% [ 0|1 (default 1)] perform breadth first search for exploration in the parameter space 

options.plcp.dfs = false;
% [ 0|1 (default 0) ] perform depth first search for exploration in the parameter space

options.plcp.debug = 0;
% [0|1|2 (default 0)] debugging level, 1-no plots, 2-including plots

options.plcp.fixedstep = false; % [ 0| 1 (default 0)]
% always perform fixed step over the facet to detect neighbors

options.plcp.maxlayers = inf;
% This option limits the number of layers to be explored only for BFS type
% of exploration in the parameter space. This is particularly of interest
% when there are too many regions far from the initial region that can be
% discarded. By default it is not limited.

options.plcp.maxregions = inf;
% Maximum number of regions to be generated. By default it is not limited.

options.plcp.QRfactor = false; %[ false by default ]
% use recursive QR factorization for pivoting (faster but numerically bad)

options.plcp.checkoverlaps = false; %[ false by default ]
% check for overlaps while exploring (significantly slows down the
% computation)

options.plcp.rescue = false; %[ false by default ]
% if the variable step approach fails to find a neighbor, retry with a
% fixed step

options.plcp.maxsteps = 200; %[ 200 ], minimum number of steps are 2
% maximum number of steps to compute with the fixed-step approach
%options.plcp.stepsize = 1e-5; % [1e-3 < stepsize <abs_tol (1e-5 by default)],
% step size to be performed with the fixed step approach

options.plcp.maxpivots = 100; %[ default 100, maximum 500 unless default recursion limit is not changed in Matlab settings ],
% maximum number of pivots to perform for searching a neighbor

options.plcp.adjcheck = 0; % [default 0, otherwise 1]
% force checking of the adjacency list inside the PLCP solver

%% ENUM_PLCP options
options.enum_plcp.maxLPs = Inf;
% maximum number of LPs to solve

options.enum_plcp.maxregions = Inf;
% maximumum number of regions allowed to discover

%% LCP options
options.lcp.zerotol = 1e-10;
% Less than this treshold the value is considered as zero

options.lcp.lextol = 1e-9;
% Lexicographic tolerance - a small treshold from which values are considered as equal.

options.lcp.maxpiv = 1e4; 
% Maximum number of pivots to be performed.          

options.lcp.nstepf = 50;
% If options.routine is 0, then every "nstepf" pivot steps the basis is refactorized to
% avoid numerical problems for LUMOD. For other routines the factorization
% is performed at each step.

options.lcp.clock = 0;
% Show the information about the computational time.                                                                          

options.lcp.verbose = 0;
% Verbose output. Show progress of pivoting algorithm including entering,
% leaving variables, actual basis and basis solution.                   

options.lcp.routine = 1;
% Routine which should be used to obtain a basis solution.                                          
% 0 - Corresponds to LUmod package that performs factorization in the form
%    LA = U. Depending on the change in A  factors L, U  are updated. This is
%    the fastest method. 
% 1 - Corresponds to DGESV simple driver from LAPACK package which solves
%    the system AX = B by factorizing A and overwriting B with the solution X.
%    Since the factorization is performed at each pivot step, this method
%    tends to be much slower than method 0.
% 2 - Corresponds to DGELS simple driver which solves overdetermined or
%    underdetermined real linear systems min ||b - Ax||_2 involving an 
%    M-by- N  matrix A, or its transpose, using a QR or LQ factorization of
%    A. Since the factorization is performed at each pivot step, this
%    method tends to be much slower than method 0.                                              

options.lcp.timelimit= 3600;
% Time limit in seconds. If this limit is exceeded,  the pivoting algorithm
% is terminated and current basis is returned.                                 

options.lcp.normalize = 1;
% Input matrices M, q  get scaled by D_1  and D_2 when invoking this
% option: M_n = D_1MD_2, q_n = D_1q, and solution is recovered as z =
% D_2z_n, w = Mz+q.                                               

options.lcp.normalizethres = 1e6;
% If the normalize option is on, then the matrix scaling is performed 
% only if 1 norm of matrix M (maximum absolute column sum) is above this
% threshold. 


%% GLPK options
options.glpk.msglev = 0; %  Level of messages output by solver routines:
% 0 - No output.
% 1 - Error messages only.
% 2 - Normal output.
% 3 - Full output (includes informational messages).

options.glpk.lpsolver = 1;
% 1 - Revised simplex method (default)
% 2 - Interior point method 
% 3 - Simplex method with exact arithmetic

options.glpk.scale = 1; 
% 0 - No scaling.
% 1 - Equilibration scaling (default)
% 2 - Geometric mean scaling, then equilibration scaling
% 3 - Geometric then Equilibrium scaling 
% 4 - Round to nearest power of 2 scaling

options.glpk.dual = 1;  
% 0 - Do not use the dual simplex
% 1 - If initial basic solution is dual feasible, use the dual simplex
% 2 - Use two phase dual simplex, or if primal simplex if dual fails

options.glpk.price = 1; % Pricing option (for both primal and dual simplex):
% 0 - Textbook pricing
% 1 - Steepest edge pricing

options.glpk.r_test = 1; % Ratio test Technique
% (default: 1). 
% 0 - stardard (textbook)
% 1 - Harris's two-pass ratio test

options.glpk.relax = 0.07;
% (default: 0.07). Relaxation parameter used in the ratio test. If it is
% zero, the textbook ratio test is used. If it is non-zero (should be
% positive), Harris two-pass ratio test is used. In the latter case on the
% first pass of the ratio test basic variables (in the case of primal
% simplex) or reduced costs of non-basic variables (in the case of dual
% simplex) are allowed to slightly violate their bounds, but not more than
% relax*tolbnd or relax*toldj  (thus, relax is a percentage of tolbnd or toldj).

options.glpk.tolbnd = 1e-7;
% (default: 1e-7). Relative tolerance used to check ifthe current basic
% solution is primal feasible. It is not recommended that you change this
% parameter unless you have a detailed understanding of its purpose. 

options.glpk.toldj = 1e-7;
% (default: 1e-7). Absolute tolerance used to check if the current basic
% solution is dual feasible. It is not recommended that you change this
% parameter unless you have a detailed understanding of its purpose. 

options.glpk.tolpiv = 1e-9;
% (default: 1e-9). Relative tolerance used to choose eligible pivotal
% elements of the simplex table. It is not recommended that you change this
% parameter unless you have a detailed understanding of its purpose.

options.glpk.round = 0; % Solution rounding option
% 0 - Report all primal and dual values "as is" (default).
% 1 - Replace tiny primal and dual values by exact zero.

options.glpk.objll = -1e12;
% (default: -DBL_MAX). Lower limit of the
% objective function. If on the phase II the objective
% function reaches this limit and continues decreasing, the
% solver stops the search. This parameter is used in the
% dual simplex method only.

options.glpk.objul = 1e12;
% (default: +DBL_MAX). Upper limit of the
% objective function. If on the phase II the objective
% function reaches this limit and continues increasing,
% the solver stops the search. This parameter is used in
% the dual simplex only.

options.glpk.itlim = 1e4;
% (default: -1). Simplex iterations limit.
% If this value is positive, it is decreased by one each
% time when one simplex iteration has been performed, and
% reaching zero value signals the solver to stop the search.
% Negative value means no iterations limit.

options.glpk.itcnt = 200;
% (default: 200). Output frequency, in iterations.
% This parameter specifies how frequently the solver sends
% information about the solution to the standard output.

options.glpk.usecuts = 1;
% (default: 1). ( for MIP only )
% glp_intopt generates and adds cutting planes to
% the MIP problem in order to improve its LP relaxation
% before applying the branch&bound method
% 0 - all cuts off
% 1 - Gomoy's mixed integer cuts
% 2 - Mixed integer rounding cuts
% 3 - Mixed cover cuts
% 4 - Clique cuts
% 5 - all cuts

options.glpk.pprocess = 2; % Pre-processing technique option ( for MIP only )
% (default: 2) 
% 0 - disable preprocessing
% 1 - perform preprocessing for root only
% 2 - perform preprocessing for all levels

options.glpk.binarize = 0; % Binarizeation option ( for mip only ):
% (default: 0 ) 
% ( used only if presolver is enabled )
% 0 - do not use binarization
% 1 - replace general integer variables by binary ones


options.glpk.tmlim = -1;
% (default: -1.0). Searching time limit, in
% seconds. If this value is positive, it is decreased each
% time when one simplex iteration has been performed by the
% amount of time spent for the iteration, and reaching zero
% value signals the solver to stop the search. Negative
% value means no time limit.

options.glpk.branch = 2; % Branching heuristic option (for MIP only)
% 0 - Branch on the first variable
% 1 - Branch on the last variable
% 2 - Branch on the most fractional variable (default)
% 3 - Branch using a heuristic by Driebeck and Tomlin

options.glpk.btrack = 2; % Backtracking heuristic option (for MIP only)
% 0 - Depth first search.
% 1 - Breadth first search.
% 2 - best local bound (default)
% 3 - Backtrack using the best projection heuristic.

options.glpk.tolint = 1e-6;
% (default: 1e-5). Relative tolerance used
% to check if the current basic solution is integer
% feasible. It is not recommended that you change this
% parameter unless you have a detailed understanding of
% its purpose.

options.glpk.outdly = 0;
% (default: 0.0). Output delay, in seconds.
% This parameter specifies how long the solver should
% delay sending information about the solution to the standard
% output. Non-positive value means no delay.

options.glpk.tolobj = 1e-7;
% (default: 1e-7). Relative tolerance used
% to check if the value of the objective function is not
% better than in the best known integer feasible solution.
% It is not recommended that you change this parameter
% unless you have a detailed understanding of its purpose.

options.glpk.presol = 1;
% (default: 1). If this flag is set, the routine
% lpx_simplex solves the problem using the built-in LP presolver.
% Otherwise the LP presolver is not used.

options.glpk.save = 0;
% (default: 0). If this parameter is nonzero save a copy of
% the original problem to file.

options.glpk.mipgap = 0; % The relative mip gap tolerance. 
% (default: 0.0)  
% If the relative mip gap for currently known best integer feasible
% solution falls below this tolerance, the solver terminates
% the search.  This allows obtaining suboptimal interger
% feasible solutions if solving the problem to optimality
% takes too long.


%% GUROBI options
% details in http://www.gurobi.com/doc/40/refman/node572.html
% only some important settings are here, the rest is kept to be default

% Determines whether dual reductions are performed in presolve. You should
% disable these reductions if you received an optimization status of
% INF_OR_UNBD and would like a more definitive conclusion.
% http://www.gurobi.com/documentation/5.6/reference-manual/dualreductions
%
% In short, setting DualReduction=0 prevents Gurobi from returning the
% INF_OR_UNBD status upon which mpt_call_gurobi() needs to solve the
% problem a second time to obtain a definite answer.
%
% The default in Gurobi 5.6.3 is DualReductions=1
options.gurobi.DualReductions = 1;

options.gurobi.BarIterLimit=Inf;
% >0, Limits the number of barrier iterations performed (barrier only).

options.gurobi.CutOff = Inf;
% If the objective value for the optimal solution is better than the
% specified cutoff, the solver will return the optimal solution. Otherwise,
% it will terminate with a CUTOFF status (see the Status Code section for
% further details). Infinity for minimization, -Infinity for maximization   

options.gurobi.IterationLimit=1e6; 
% >0, Limits the number of simplex iterations performed.

options.gurobi.NodeLimit=Inf;
% >0, Limits the number of MIP nodes explored (MIP only).

options.gurobi.SolutionLimit=Inf;
% >0, Limits the number of feasible solutions found (MIP only).

options.gurobi.TimeLimit=Inf;
% >0, Limits the total time expended (in seconds).

options.gurobi.BarConvTol = 1e-8;
% [1e-10, 1] Barrier convergence tolerance (barrier only). The barrier solver
% terminates when the relative difference between the primal and dual
% objective values is less than the specified tolerance.  

options.gurobi.BarQCPConvTol=1e-6;
% [0, 1] The barrier solver terminates when the relative difference between the
% primal and dual objective values is less than the specified tolerance
% (with a GRB_OPTIMAL status). Tightening this tolerance may lead to a more
% accurate solution, but it may also lead to a failure to converge.   

options.gurobi.FeasibilityTol=1e-6;
% [1e-9,1e-2] Primal feasibility tolerance. All constraints must be
% satisfied to a tolerance of FeasibilityTol. 

options.gurobi.IntFeasTol=1e-5;
% [1e-9, 1e-1] Integer feasibility tolerance (MIP only). An integrality restriction on a
% variable is considered satisfied when the variable's value is less than
% INTFEASTOL from the nearest integer value.  

options.gurobi.OptimalityTol=1e-6;
% [1e-9, 1e-2], Dual feasibility tolerance. Reduced costs must all be smaller than
% OptimalityTol in the improving direction in order for a model to be
% declared optimal.  

options.gurobi.MIPGap=1e-4;
% >0, Relative MIP optimality gap (MIP only). The MIP engine will terminate
% (with an optimal result) when the gap between the lower and upper
% objective bound is less than MIPGap times the upper bound.  

options.gurobi.PSDTol = 1e-6;
% >0, Positive semi-definite tolerance (QP/MIQP only). Sets a limit on the
% amount of diagonal perturbation that the optimizer is allowed to perform
% on the Q matrix in order to correct minor PSD violations. If a larger
% perturbation is required, the optimizer will terminate with an
% GRB_ERROR_Q_NOT_PSD error.     

options.gurobi.InfUnbdInfo=0;
%[0 or 1] Determines whether simplex (and crossover) will compute
%additional information when a model is determined to be infeasible or
%unbounded. Set this parameter if you want to query the unbounded ray for
%unbounded models (through the UnbdRay attribute), or the infeasibility
%proof for infeasible models (through the FarkasDual and FarkasProof
%attributes).      

options.gurobi.NormAdjust = -1;
%Chooses from among multiple pricing norm variants. The details of how this
%parameter affects the simplex pricing algorithm are subtle and difficult
%to describe, so we've simply labeled the options 0 through 3. The default
%value of -1 chooses automatically.    

options.gurobi.PerturbValue = 0.0002;
% [0,0.01] Default 0.0002. Magnitude of the simplex perturbation. Note that
% perturbation is only applied when progress has stalled, so the parameter
% will often have no effect.   

options.gurobi.ScaleFlag = 1;
%[0, 1] Default 1.  Enables or disables model scaling. Scaling usually
%improves the numerical properties of the model, which typically leads to
%reduced solution times, but it may sometimes lead to larger constraint
%violations in the original, unscaled model.   

options.gurobi.ObjScale=0;
% >-1,  Divides the model objective by the specified value to avoid
% numerical errors that may result from very large objective coefficients.
% The default value of 0 decides on the scaling automatically. A value less
% than zero uses the maximum coefficient to the specified power as the
% scaling (so ObjScale=-0.5 would scale by the square root of the largest
% objective coefficient).     

options.gurobi.BarCorrectors = -1;
%>-1. Default -1. Limits the number of central corrections performed in each barrier iteration

options.gurobi.Method=1;
% Algorithm used to solve continuous models or the root node of a MIP model
% (-1=automatic, 0=primal simplex, 1=dual simplex, 2=barrier, 3=concurrent,
% 4=deterministic concurrent).Simplex algorithm (0=primal, 1=dual). 
    
options.gurobi.Presolve=-1;
% Controls the presolve level. A value of -1 corresponds to an automatic
% setting. Other options are off (0), conservative (1), or aggressive (2).

options.gurobi.TimeLimit=Inf;
% >0, Limits the total time expended (in seconds).

options.gurobi.Threads=0;
% [0, number of processors]. Controls the number of threads to apply to parallel MIP. The default
% value of 0 sets the thread count equal to the maximum value, which is the
% number of processors in the machine.  

options.gurobi.OutputFlag = 0;
% Enables or disables solver output. Use LogFile and LogToConsole for
% finer-grain control. Setting OutputFlag to 0 is equivalent to setting
% LogFile to "" and LogToConsole to 0.   

options.gurobi.DisplayInterval=5;
% (default = 5). Controls the frequency at which log lines are printed (in seconds).

options.gurobi.Aggregate = 0;
% Enables or disables aggregation in presolve. In rare instances,
% aggregation can lead to an accumulation of numerical errors. Turning it
% off can sometimes improve solution accuracy.

%% MOSEK options
try
    [~, res] = mosekopt('param');
    options.mosek = res.param;
end

%% NAG options
options.nag.qp.ftol = 1e-9;
% (default 1e-9)
% the maximum acceptable violation in each constraint at a "feasible" point.

options.nag.qp.rank_tol = 1e-20;
% (default = 1e-13)
% enables the user to control the condition number of the triangular factor R

options.nag.qp.crash_tol = 0.1; 
% (default = 0.01)
% a constraint of the form a'*x >= l will be included in the initial
% working set if |a'*x-l| <= crash_tol times (1+|l|) 0.0 <
% options.crash_tol <= 1.0 options.nag.orthog = 0;
% DO NOT PROVIDE "crash_tol" as 1 because it causes NAG to segfault

options.nag.qp.reset_ftol = 5;
% (default = 5)
% this option is part of an anti-cycling procedure designed to guarantee
% progress even on highly degenerate problems 

options.nag.qp.max_iter = 1e6;
% maximum number of iterations to be performed

options.nag.qp.fcheck = 50;
% every fcheck iterations, a numerical test is
% made to see if the current solution x satisfies
% the constraints in the working set

options.nag.qp.inf_bound = 1e12;
% defines the "infinite" bound in the definition
% of the problem constraints

options.nag.lp.ftol = 1e-9;
% the maximum acceptable violation in each constraint at a "feasible" point.

options.nag.lp.optim_tol = 1e-13;
% enables the user to control the condition number of the triangular factor R

options.nag.lp.crash_tol = 0.01; 
% a constraint of the form a'*x >= l will be included in the initial
% working set if |a'*x-l| <= crash_tol times (1+|l|) 0.0 <
% options.crash_tol <= 1.0 options.nag.orthog = 0;

options.nag.lp.reset_ftol = 5;
% this option is part of an anti-cycling procedure designed to guarantee
% progress even on highly degenerate problems 

options.nag.lp.max_iter = 1e6;
% maximum number of iterations to be performed

options.nag.lp.fcheck = 50;
% every fcheck iterations, a numerical test is
% made to see if the current solution x satisfies
% the constraints in the working set

options.nag.lp.inf_bound = 1e12;
% defines the "infinite" bound in the definition
% of the problem constraints


%% QPC
options.qpip.mu = 0.0; 
% (default mu=0.0). Desired complementarity gap target (point on central-path)

options.qpip.method = 0; 
% If method=1 (default), then a faster but less accurate linear
% solve step is used. Conversely, if method=0 then a slower but
% more accurate linear solve step is used.

%% QPspline
options.qpspline.maxiter = 10000;
% maximum iterations

options.qpspline.abs_tol = 1e-8; 
% absolute tolerance

options.qpspline.nstepf = 30;
% after these n steps the basis will be refactored with Q,R, factors

options.qpspline.nqrelems = 20;
% Do recursive QR factorization if the number of changed elements
% (rows/cols) is less than this treshold. If the treshold is large, the
% recursive factorization may actually consume more time than direct
% factorization. Large treshold causes also accumulation of numerical
% errors. 

options.qpspline.timelimit = 3600;
% time limit on the computations in seconds

options.qpspline.verbose = 0;
% show the iteration progress

%% QUADPROG
if exist('mskoptimset', 'file')
    % MOSEK
    % I really hate doing this. Hey, MOSEK, stop messing with quadprog!
    options.quadprog = mskoptimset('quadprog');
elseif exist('quadprog','file')==2
    options.quadprog = optimset('quadprog');
end
options.quadprog.Display = 'off';
% no display

options.quadprog.MaxIter = 1e6;
% Maximum number of iterations allowed

options.quadprog.TolFun = 1e-10;
% Termination tolerance on the function value

options.quadprog.TolX = 1e-10;
% Termination tolerance on X

options.quadprog.Algorithm = 'trust-region-reflective';
%options.quadprog.Algorithm = 'interior-point-convex';
% Algorithm to use

options.quadprog.LargeScale='off';
% do not use largescale algorithm

%% LINPROG
if exist('mskoptimset', 'file')
    % MOSEK
    % I really hate doing this. Hey, MOSEK, stop messing with linprog!
    options.linprog = mskoptimset('linprog');
elseif exist('linprog','file')==2
    options.linprog = optimset('linprog');
end

options.linprog.Display = 'off';
% no display

options.linprog.MaxIter = 1e6;
% Maximum number of iterations allowed

options.linprog.TolFun = 1e-10;
% Termination tolerance on the function value

options.linprog.TolX = 1e-10;
% Termination tolerance on X

%% SEDUMI
options.sedumi.fid = 0; 
% silent output

options.sedumi.alg = 2;
% If alg=0, then a first-order wide
% region algorithm is used, not recommended. If alg=1, then SeDuMi uses
% the centering-predictor-corrector algorithm with v-linearization.
% If alg=2, then xz-linearization is used in the corrector, similar
% to Mehrotra's algorithm. The wide-region centering-predictor-corrector
% algorithm was proposed in Chapter 7 of
% J.F. Sturm, Primal-Dual Interior Point Approach to Semidefinite Pro-
% gramming, TIR 156, Thesis Publishers Amsterdam, 1997.

options.sedumi.theta = 0.85;
options.sedumi.beta = 0.5;
% By default, theta=0.25 and beta=0.5. These
% are the wide region and neighborhood parameters. Valid choices are
% 0 < theta <= 1 and 0 < beta < 1. Setting theta=1 restricts the iterates
% to follow the central path in an N_2(beta)-neighborhood.

options.sedumi.stepdif = 2;
options.sedumi.w = [1 1];
% By default, stepdif = 2 and w=[1 1].
% This implements an adaptive heuristic to control step-differentiation.
% You can enable primal/dual step length differentiation by setting stepdif=1 or 0.
% If so, it weights the rel. primal, dual and gap residuals as
% w(1):w(2):1 in order to find the optimal step differentiation.

options.sedumi.eps = 1e-10;
% The desired accuracy. Setting pars.eps=0 lets SeDuMi run
% as long as it can make progress. By default eps=1e-8.

options.sedumi.bigeps = 1e-5;
% In case the desired accuracy pars.eps cannot be achieved,
% the solution is tagged as info.numerr=1 if it is accurate to pars.bigeps,
% otherwise it yields info.numerr=2.

options.sedumi.maxiter = 1e6;
% maximum number of iterations

%  Various parameters for controling the Preconditioned conjugate
%  gradient method (CG), which is only used if results from Cholesky are inaccurate.
options.sedumi.cg.maxiter = 49;
% (default 49)
% Maximum number of CG-iterates (per solve). Theoretically needed
% is |add|+2*|skip|, the number of added and skipped pivots in Cholesky.

options.sedumi.cg.restol = 5e-3;
% (default 5E-3)
% Terminates if residual is a "cg.restol" fraction of duality gap.
% Should be smaller than 1 in order to make progress.

options.sedumi.cg.refine = 1;
% (default 1)
% Number of refinement loops that are allowed. The maximum number
% of actual CG-steps will thus be 1+(1+cg.refine)*cg.maxiter. 

options.sedumi.cg.stagtol = 4e-14;
% (default 5E-14)
% Terminates if relative function progress less than stagtol.

options.sedumi.cg.qprec = 0;
% (default 0)
% Stores cg-iterates in quadruple precision if qprec=1.

% Various parameters for controling the Cholesky solve.
options.sedumi.chol.canceltol = 1e-12;
% (default 1E-12)
% Rel. tolerance for detecting cancelation during Cholesky 

options.sedumi.chol.maxu = 5e5;
% (default 5E5)
% Adds to diagonal if max(abs(L(:,j))) > chol.maxu otherwise.

options.sedumi.chol.abstol = 1e-20;
% (default 1e-20)
% Skips pivots falling below abstol.

options.sedumi.chol.maxuden = 5e2;
% (default 5E2)
% pivots in dense-column factorization so that these factors
% satisfy max(abs(Lk)) <= maxuden .

options.sedumi.errors = 0;
% (default 0)
% If this field is 1 then SeDuMi outputs some error
% measures as defined in the Seventh DIMACS Challenge. 

%% mpt_enumpqp
%     .prune_infeasible: if true, list of candidate active sets is pruned
%                        by removing entries that were previously
%                        discovered as infeasible (default=true)
options.enum_pqp.prune_infeasible = true;

%         .feasible_set: if 'regions' (default), the feasible set is
%                        identical to critical regions. If 'projection', 
%                        the set of feasible parameters is computed by
%                        projection. If 'outerbox', the feasible set is
%                        costructed as the outer box approximation of the
%                        union of critical regions. If 'rn', the feasible
%                        set is returned as R^n. (default='regions')
options.enum_pqp.feasible_set = 'regions';

% if .remove_redundant=true, redundant inequalities will be removed from
% the pQP problem (improves runtime, but can be slow for large problems)
options.enum_pqp.remove_redundant = true;

% If .accept_regions='fulldim', only fully dimensional regions are returned.
% If .accept_regions='nonempty', also flat (lower-dimensiona) regions will
% be returned
options.enum_pqp.accept_regions = 'fulldim';

% Backoff for inequality constraints in IPDPolyhedron
options.enum_pqp.ineq_backoff = 1e-9;
