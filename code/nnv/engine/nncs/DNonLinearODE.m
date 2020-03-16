classdef DNonLinearODE < handle
    % Discrete time nonlinear ODE class. Grapper of nonlinearSysDT in CORA.
    % Chelsea Sidrane April 25, 2019 && Dung Tran: 4/29/2019
    % Based on NonLinearODE.m
    
    properties
        options = []; % option for recahable set computation
        dynamics_func = []; % function to describe dynamics of the system
        dim = 0; % dimension of the system (number of state variable)
        nI = 0; % number of control inputs
        nO = 0; % number of outputs
        C; % output matrix y = Cx
        Ts = 0 ; % sampling time
        intermediate_reachSet = []; % intermediate reachable set between steps
        % if discrete timestep of plant is faster than discrete timestep of
        % controller
    end
    
    methods
        % constructor
        function obj = DNonLinearODE(dim, nI, dynamics_func, Ts, outputMat)
            % construct DnonLinearODE plant
            % @dim: dimension of the plant
            % @nI: number of input
            % @dynamics_func: dynamics of the plant, the input should have
            % @ character to specify function, for example @car_dynamics
            % Ts: timestep
            % @outputMat: output matrix
            
            % Note: we construct a DnonLinearODE plant with default option
            
            if dim < 1
                error('Dimension should be >=1');
            end
            
            if nI < 0
                error('Number of inputs should be >= 0');
            end
            
            obj.dim = dim;
            obj.nI = nI;
            obj.dynamics_func = dynamics_func;
            obj.Ts = Ts;
                        
            % default option
            option.tStart = 0;
            option.tFinal = Ts; % we provide method to change the option
            
            option.x0 = [];
            option.R0 = []; % initial set for reachability analysis
            option.timeStep = Ts; % time step for reachable set computation
            option.taylorTerms = 4; % number of taylor terms for reachable sets
            option.zonotopeOrder = 2; % zonotope order
            option.intermediateOrder = 5; 
            option.reductionTechnique = 'girard';
            option.errorOrder = 1;
            option.polytopeOrder = 2; % polytope order
            option.reductionInterval = 1e3;
            option.maxError = 0.5*ones(dim, 1);
            option.advancedLinErrorComp = 0;
            option.tensorOrder=2;
            option.uTrans = 0;
            
            obj.options = option; % default option
            obj.set_output_mat(outputMat);
            
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
        
        % set polytope Order
        function set_polytopeOrder(obj, polytopeOrder)
            
            if polytopeOrder < 1
                error('Invalid polytope order');
            end
            
            obj.options.polytopeOrder = polytopeOrder;
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
        
        % set originContained
        function set_originContained(obj, originContained)
            obj.options.originContained = originContained;
        end
        
        % set advancedLinErrorComp
        function set_advancedLinErrorComp(obj, advancedLinErrorComp)
            obj.options.advancedLinErrorComp = advancedLinErrorComp;
        end
        
        % set plot Type
        function set_plotType(obj, plotType)
            obj.options.plotType = plotType;
        end
        
        % set uTrans
        function set_uTrans(obj, uTrans)
            obj.options.uTrans = uTrans;
        end
        
        % set input set U
        function set_U(obj, U)
            obj.options.U = U;
        end
        
        % set initial set
        function set_R0(obj, R0)
            obj.options.R0 = R0;
        end
        
        % set tFinal
        function set_tFinal(obj, tFinal)
            if tFinal <= 0
                error('tFinal should be > 0');
            end
            obj.options.tFinal = tFinal;
        end
        
        % set tStart
        function set_tStart(obj, tStart)
            
            if tStart < 0
                error('Invalid tStart');
            end
            
            obj.options.tStart = tStart;
        end
        
        % set time step 
        function set_timeStep(obj, timeStep)
            if timeStep <= 0
                error('Invalid time step');
            end
            obj.options.timeStep = timeStep;
        end
        
        % set inital state for simulation x0
        function set_x0(obj, x0)
            if length(x0) ~= obj.dim
                error('Dimension mismatch between initial state and the system');
            end
            
            obj.options.x0 = x0;
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
    
    % reachability anlaysis methods
    methods
        
         % reachability analysis using zonotope
        function [R, reachTime] = reach_zono(obj, init_set, input_set, tFinal)
           % @init_set: initial set of state
           % @input_set: input set u
           % @timeStep: time step in reachable set computation
           
           % this is a grapper of reach method for DT nonlinear system in CORA
           % the initial set and input set are Zono in nnv
           
           start = tic;
           if ~isa(init_set, 'Zono')
               error('Initial set is not a zonotope');
           end
           
           if ~isa(input_set, 'Zono')
               error('Input set is not a zonotope');
           end
           
           % zonotope is a CORA class
           R0 = zonotope([init_set.c, init_set.V]);
           U = zonotope([input_set.c, input_set.V]);
           
           obj.set_R0(R0);
           obj.set_U(U);
                                
           sys = nonlinearSysDT(obj.dim, obj.nI, obj.dynamics_func, obj.options); % CORA nonlinearSys class
           R = reach(sys, obj.options); % CORA reach method using zonotope and conservative linearization
                     
           reachTime = toc(start);
                   
        end
        
       % step reach using star set used for neural network control system
        function S = stepReachStar(obj, init_set, input_set)
            % @init_set: initial set, a star
            % @input_set: input set, a star
            % @R: reachable set at the timeStep, a star
            
            if ~isa(init_set, 'Star')
                error('Initial set is not a star');
            end
            
            if ~isa(input_set, 'Star')
                error('Input set is not a star');
            end
            
            I = init_set.getZono;
            U = input_set.getZono;
            
            [R, ~] = obj.reach_zono(I, U, obj.options.tFinal);
            
            N = length(R);           
            Z = R{N,1}; % get the last zonotope in the reach set
            Z = Z.Z; % get c and V 
            c = Z(:,1); % center vector
            V = Z(:, 2:size(Z, 2)); % generators
            
            % Zono is the verivital zonotope class
            Z = Zono(c, V);
            S = Z.toStar;
            
            for i=1:N
                Z = R{i,1}; 
                Z = Z.Z; % get c and V 
                c = Z(:,1); % center vector
                V = Z(:, 2:size(Z, 2)); % generators

                Z = Zono(c, V);
                S = Z.toStar;
                obj.intermediate_reachSet = [obj.intermediate_reachSet S];
            end
            
            % the last zonotope in the reach set is returned
            
        end
        
        % evaluate (simulate) the plant with specific input and state
        % using ode45 solver
        function y = evaluate(obj, x0, u)
            % @x0: initial state
            % @u: control input
            % @tspan: time points
            
            % author: Dung Tran
            % date: 1/29/2019
            
            y = obj.dynamics_func([], x0, u, []);          
            
        end
        
        % implement box?
        
        % implement polyhedron?
        
        % add simulation based evaluation?
        
    end
end