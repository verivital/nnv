classdef NonLinearODE < handle
    % Nonlinear ODE class is a grapper of nonlinearSys class in CORA
    %   Dung Tran: 11/19/2018
    % Last Revision: September 3, 2020 - Diego Manzanas
    %      - CORA 2020 updates
    %
    %
    % obj = NonLinearODE(dim, nI, dynamics_func, reachTimeStep, controlPeriod, outputMat)
            % obj = construct nonLinearODE plant
            % @dim: dimension of the plant
            % @nI: number of input
            % @dynamics_func: dynamics of the plant, the input should have
            % @ character to specify function, for example @car_dynamics
            % @reachTimeStep: reachability step for the plant
            % @controlPeriod: control period for the plant
            % @outputMat: output matrix
    
    properties
        options = []; % option for recahable set computation
        params = [];
        dynamics_func = []; % function to describe dynamics of the system
        dim = 0; % dimension of the system (number of state variable)
        nI = 0; % number of control inputs
        nO = 0; % number of outputs
        C; % output matrix y = Cx
        controlPeriod = 0.1; % control period
        intermediate_reachSet = []; % intermediate reachable set between steps 
        intermediate_pointSet = []; % intermediate reachable set between steps at a single point in time
        reachStep = 0.01; % reachability step for the plant
        cora_set = []; % avoid overapproximation from polyzono to zono, use cora set from last reach set
    end
    
    methods
        
        % constructor from a matlab function reference in dynamics_func
        function obj = NonLinearODE(dim, nI, dynamics_func, reachTimeStep, controlPeriod, outputMat)
            % construct nonLinearODE plant
            % @dim: dimension of the plant
            % @nI: number of input
            % @dynamics_func: dynamics of the plant, the input should have
            % @ character to specify function, for example @car_dynamics
            % @reachTimeStep: reachability step for the plant
            % @controlPeriod: control period for the plant
            % @outputMat: output matrix
            
            % author: Dung Tran
            % date: 11/19/2018
            % update: 3/15/2020
            
            % Note: we construct a nonLinearODE plant with default option

            if dim < 1
                error('Dimension should be >=1');
            end
            
            if nI < 0
                error('Number of inputs should be >= 0');
            end
            
            obj.dim = dim;
            obj.nI = nI;
            obj.dynamics_func = dynamics_func;
            
            obj.NonLinearODE_init();
            obj.set_timeStep(reachTimeStep);
            obj.reachStep = reachTimeStep;
            obj.set_tFinal(controlPeriod);
            obj.controlPeriod = controlPeriod;
            obj.set_output_mat(outputMat);
        end
        
        function obj = NonLinearODE_hyst(hyst_ha)
            % pull arguments from hyst ha
            obj = NonLinearODE();
            obj.ha = hyst_ha;
        end
        
        % set default parameters
        function NonLinearODE_init(obj)
            % default option
            param.tStart = 0;
            param.tFinal = 0.1; % we provide method to change the option
            
            param.R0 = []; % initial set for reachability analysis
            option.timeStep = 0.01; % time step for reachable set computation
            option.taylorTerms = 4; % number of taylor terms for reachable sets
            option.zonotopeOrder = 20; % zonotope order
            option.intermediateOrder = 20; 
            % option.reductionTechnique = 'girard';
            option.errorOrder = 1;
            option.reductionInterval = 1e3;
            option.maxError = 0.1*ones(obj.dim, 1);
            option.tensorOrder=3; % Recommended 2 or 3 (minimum 3 for poly)
            option.alg = 'lin'; % 'lin-adaptive' or 'poly-adaptive' recommended (no need to select other parameters)
            
            obj.options = option; % default option
            obj.params = param; % default parameters
        end
        
        % set taylor terms
        function set_taylorTerms(obj, taylorTerms)
            if taylorTerms < 1
                error('Invalid Taylor Terms');
            end
            obj.options.taylorTerms = taylorTerms;
        end
        
        % set zonotope order
        function set_zonotopeOrder(obj, zonotopeOrder)
            if zonotopeOrder < 1
                error('Invalid zonotope order');
            end
            obj.options.zonotopeOrder = zonotopeOrder;
        end
        
        % set intermediate order
        function set_intermediateOrder(obj, intermediateOrder)
            if intermediateOrder < 1
                error('Invalid intermediate order');
            end
            obj.options.intermediateOrder = intermediateOrder;
        end
        
        % set reduction technique
        function set_reductionTechnique(obj, reductionTechnique)
            obj.options.reductionTechnique = reductionTechnique;
        end
        
        % set error order
        function set_errorOrder(obj, errorOrder)
            if errorOrder < 1
                error('Invalid error order');
            end
            obj.options.errorOrder = errorOrder;
        end
        
        % set reduction Interval
        function set_reductionInterval(obj, reductionInterval)
            if reductionInterval < 0
                error('Invalid reduction interval');
            end
            obj.options.reductionInterval = reductionInterval;
        end
        
        % set max Error       
        function set_maxError(obj, maxError)
            if size(maxError, 1) ~= obj.dim && size(maxError, 2) ~= 1
                error('Invalid max error');
            end
            obj.options.maxError = maxError;
        end
        
        % set tensor Order
        function set_tensorOrder(obj, tensorOrder)
            if tensorOrder < 1
                error('Invalid tensor Order');
            end
            obj.options.tensorOrder = tensorOrder;
        end
        
        % set input set U
        function set_U(obj, U)
            obj.params.U = U;
        end
        
        % set initial set
        function set_R0(obj, R0)
            obj.params.R0 = R0;
        end
        
        % set tFinal
        function set_tFinal(obj, tFinal)
            if tFinal <= 0
                error('tFinal should be > 0');
            end
            obj.params.tFinal = tFinal;
        end
        
        % set tStart
        function set_tStart(obj, tStart)
            if tStart < 0
                error('Invalid tStart');
            end
            obj.params.tStart = tStart;
        end
        
        % set time step 
        function set_timeStep(obj, timeStep)
            if timeStep <= 0
                error('Invalid time step');
            end
            obj.options.timeStep = timeStep;
        end
        
        % set output matrix
        function set_output_mat(obj, output_mat)
            if size(output_mat, 2) ~= obj.dim
                error('Invalid output matrix');
            end
            obj.nO = size(output_mat, 1);
            obj.C = output_mat;
        end
    end
    
    
    % reachability anlaysis method
    methods
        
        % reachability analysis using zonotope
        function [R, reachTime] = reach_zono(obj, init_set, input_set, timeStep, tFinal)
           % @init_set: initial set of state
           % @input_set: input set u
           % @timeStep: time step in reachable set computation
           % @num_steps: number of steps in reachable set computation
           
           % this is a grapper of reach method for nonlinear system in CORA
           % the initial set and input set are Zono in nnv
           
           % author: Dung Tran
           % date: 11/20/2018
           
           start = tic;
           if ~isa(init_set, 'Zono')
               error('Initial set is not a zonotope');
           end
           
           if ~isa(input_set, 'Zono')
               error('Input set is not a zonotope');
           end
           
           R0 = zonotope([init_set.c, init_set.V]);
           U = zonotope([input_set.c, input_set.V]);
           
           obj.set_R0(R0);
           obj.set_U(U);
           obj.set_timeStep(timeStep);
           obj.set_tFinal(tFinal);
           
           sys = nonlinearSys(obj.dynamics_func, obj.dim, obj.nI); % CORA nonlinearSys class
           R = reach(sys, obj.params, obj.options); % CORA reach method using zonotope and conservative linearization
                     
           reachTime = toc(start);
                   
        end
        
        % reachability analysis using polynomial zonotopes
        function [R, reachTime] = reach_polyZono(obj, init_set, input_set, timeStep, tFinal)
           % @init_set: initial set of state
           % @input_set: input set u
           % @timeStep: time step in reachable set computation
           % @num_steps: number of steps in reachable set computation
           
           % this is a grapper of reach method for nonlinear system in CORA
           % the initial set and input set are Zono in nnv
           
           % author: Diego Manzanas
           % date: 06/07/2021
           
           start = tic;
           if ~isa(init_set, 'Zono')
               error('Initial set is not a zonotope');
           end
           
           if ~isa(input_set, 'Zono')
               error('Input set is not a zonotope');
           end
           
           % Set initial states
           if isempty(obj.cora_set) % set CORA params with initial state from NNV
               R0 = zonotope([init_set.c, init_set.V]);
               obj.set_R0(polyZonotope(R0));
           else % reuse last controlPeriod's reach set from CORA (reduce over-approx)
               R0 = obj.cora_set{end}.timePoint.set{end};
               obj.set_R0(R0);
           end
           % Set input set
           U = zonotope([input_set.c, input_set.V]);
           obj.set_U(U);
           obj.set_timeStep(timeStep);
           obj.set_tFinal(tFinal);
           
           % special settings for polynomial zonotopes
           polyZono.maxDepGenOrder = 50;
           polyZono.maxPolyZonoRatio = 0.01;
           polyZono.restructureTechnique = 'reducePca';
           obj.options.polyZono = polyZono;
           
           sys = nonlinearSys(obj.dynamics_func, obj.dim, obj.nI); % CORA nonlinearSys class
           R = reach(sys, obj.params, obj.options); % CORA reach method
                     
           reachTime = toc(start);
        end
        
        % step reach using star set used for neural network control system
        function S = stepReachStar(obj, init_set, input_set, varargin)
            % @init_set: initial set, a star
            % @input_set: input set, a star
            % @R: reachable set at the timeStep, a star
            
            % author: Dung Tran
            % date: 11/20/2018
            % Update: Diego Manzanas - 06/07/2021
            %     - Add 'poly' reachability analysis (varargin)
            
            % Check inputs
            if ~isa(init_set, 'Star')
                error('Initial set is not a star');
            end
            
            if ~isa(input_set, 'Star')
                error('Input set is not a star');
            end
            
            I = init_set.getZono;
            U = input_set.getZono;
            
            if ~isempty(varargin)
                if string(varargin{1}) == "poly" || string(varargin{1}) == "lin" || string(varargin{1}) == "lin-adaptive" || string(varargin{1}) == "poly-adaptive"
                    obj.options.alg = varargin{1};
                else
                    error('Incorrect reach variable options')
                end
            end
            
            % Perform reachability
            if contains(obj.options.alg,'lin')
                [R, ~] = obj.reach_zono(I, U, obj.options.timeStep, obj.params.tFinal);
            else
                [R, ~] = obj.reach_polyZono(I, U, obj.options.timeStep, obj.params.tFinal);
            end

            % Store reach sets
           if isempty(obj.cora_set)
               obj.cora_set{1} = R;
           else
               n = length(obj.cora_set);
               obj.cora_set{n+1} = R;
           end

            % Post-process CORA reach sets at t=cP
            S = [];
            cP = obj.controlPeriod;
            Rf = find((R),'time',cP);
            for mp=1:length(Rf)
                Z = Rf(mp).timePoint.set;
                Nz = length(Z); % number of sets in the computed set at cP
                for ik=1:Nz
                    try
                        Z1 = zonotope(Z{ik}); % ensure it's a zonotope
                        Z1 = Z1.Z;
                    catch
                        Z1 = zonotope(Z);
                        Z1 = Z1.Z; % get c and V 
                    end
                    c = Z1(:,1); % center vector
                    V = Z1(:, 2:size(Z1, 2)); % generators
                    Zz = Zono(c, V);
                    S = [S Zz.toStar];
                end
            end
            % we return S as the reach set at t = controlperiod
        end
        
        % evaluate (simulate) the plant with specific input and state
        function [t, y] = evaluate(obj, x0, u)
            % @x0: initial state
            % @u: control input
            % @tspan: time points
            % 
            % author: Dung Tran
            % date: 1/29/2019
            % Last Revision: 
            % 09/07/2020 - Diego Manzanas
            
            simOpt = obj.options;
            simOpt.x0 = x0;
            simOpt.tFinal = obj.params.tFinal;
            simOpt.u = u;
            odeOptions = odeset('RelTol',1e-8,'AbsTol',1e-8);
            sys = nonlinearSys(obj.dynamics_func, obj.dim, obj.nI);
            [t,y] = simulate(sys, simOpt,odeOptions);
        end

    end

    methods % other
        
        % Get interval intermediate sets from CORA (interval length of reachstep)
        function get_interval_sets(obj)
            % Get interval reach set
            Rcell = obj.cora_set;
            for r=1:length(Rcell)
                R = Rcell{r};
                N = length(R); % number of reachsets in the computation
                for i=1:N
                    Z = R(i).timeInterval.set; % get the interval reachset 
                    Nn = length(Z); % number of sets in the reachset (1 x timeStep)
                    Ss = [];
                    for ik=1:Nn
                        Z1 = Z{ik};
                        Nz = length(Z1);
                        for iz=1:Nz
                            try
                                Z2 = zonotope(Z1{iz}); % ensure it's a zonotope
                                Z2 = Z2.Z;
                            catch
                                Z2 = zonotope(Z1); % ensure it's a zonotope
                                Z2 = Z2.Z; % get c and V 
                            end
                            c = Z2(:,1); % center vector
                            V = Z2(:, 2:size(Z2, 2)); % generators
                            Zz = Zono(c, V);
                            Ss = [Ss Zz.toStar];
                        end
                    end
                    obj.intermediate_reachSet = [obj.intermediate_reachSet Ss];
                end
            end
        end
        
        % Get point intermediate sets from CORA (at each reach step)
        function get_point_sets(obj)
            % Get point reach set
            Rcell = obj.cora_set;
            for r=1:length(Rcell)
                R = Rcell{r};
                N = length(R); % number of reachsets in the computation
                for i=1:N
                    Z = R(i).timePoint.set; % get the interval reachset 
                    Nn = length(Z); % number of sets in the reachset (1 x timeStep)
                    Ss = [];
                    for ik=1:Nn
                        Z1 = Z{ik};
                        Nz = length(Z1);
                        for iz=1:Nz
                            try
                                Z2 = zonotope(Z1{iz}); % ensure it's a zonotope
                                Z2 = Z2.Z;
                            catch
                                Z2 = zonotope(Z1); % ensure it's a zonotope
                                Z2 = Z2.Z; % get c and V 
                            end
                            c = Z2(:,1); % center vector
                            V = Z2(:, 2:size(Z2, 2)); % generators
                            Zz = Zono(c, V);
                            Ss = [Ss Zz.toStar];
                        end
                    end
                    obj.intermediate_pointSet = [obj.intermediate_pointSet Ss];
                end
            end
            
        end
    
    end
    
end

