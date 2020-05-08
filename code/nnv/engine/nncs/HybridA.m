classdef HybridA < handle
    %HYBRID AUTOMATA class is a grapper of Hybrid Automaton class in CORA 
    %   Author: Diego Manzanas 05/07/2020
    %   Last revision: Diego Manzanas 05/07/2020
    % 
    %HYBRIDA Construct an instance of this class
    %      HA = HybridA(dum, nI,dynamics_func, modes);
    %      HA = NNV HA grapper class for CORA Hybrid Automaton class
    %      dim = dynamics dimensions for each mode
    %      nI = number of inputs
    %      dynamics_func = hybrid dynamics in CORA format
    %      modes = number of modes of the hybrid automata dynamics
    
    properties
        options = []; % option for recahable set computation
        dynamics_func = []; % function to describe dynamics of the system
        modes = 0; % number of HA modes
        dim = 0; % dimension of the system (number of state variable)
        nI = 0; % number of control inputs
        nO = 0; % number of outputs
        C; % output matrix y = Cx
        controlPeriod = 0.1; % control period
        intermediate_reachSet = []; % intermediate reachable set between steps
        reachStep = 0.01; % reachability step for the plant
        sysCORA;
        % this is used to store all intermediate reachable sets
    end
    
    methods
        
        function obj = HybridA(dim, nI, dynamics_func, modes)
            %HYBRIDODE Construct an instance of this class
            %   HA = HybridA(dim, nI, dynamics_func, modes);
            %   HA = NNV HA grapper class for CORA Hybrid Automaton class
            %   dim = dynamics dimensions for each mode
            %   nI = number of inputs
            %   dynamics_func = hybrid dynamics in CORA format
            %   modes = number of modes of the hybrid automata dynamics
            obj.dim = dim;
            obj.nI = nI;
%             obj.dynamics_func = @dynamics_func;
            obj.modes = modes;
            obj.C = eye(dim); % default output matrix
            obj.sysCORA = dynamics_func;
            obj.HybridA_init;
        end
        
        function HybridA_init(obj)
            % Reachability Time default option
            option.tStart = 0;
            option.tFinal = 0.1; % we provide method to change the option
            % Reachability mode location default option
            option.startLoc = 1; %initial location
            option.finalLoc = 0; %0: no final location
            
            option.x0 = []; %initial state for simulation
            option.R0 = []; %initial state for reachability analysis
            for i=1:obj.modes
                option.timeStepLoc{i} = 0.05; %time step size for reachable set computation in location i
                % Initialize inputs (no inputs)
                option.uLoc{i} = 0;
                option.uLocTrans{i} = 0;
                option.Uloc{i} = zonotope(0);
            end
            option.taylorTerms = 10;
            option.zonotopeOrder = 20;
            option.polytopeOrder = 10;
            option.errorOrder = 2;
            option.reductionTechnique = 'girard';
            option.isHyperplaneMap = 0;
            option.guardIntersect = 'polytope';
            option.enclosureEnables = 5; %choose enclosure method(s)
            option.originContained = 0;
            option.maxProjectionError = inf;
            option.linAlg = 1;
            obj.options = option; % default options
        end
        
        % set linAlg
        function set_linAlg(obj,a)
            obj.options.linAlg = a;
        end
        
        % set max projection error
        function set_maxProjectionError(obj,a)
            obj.maxProjectionError = a;
        end
        % set output matrix
        function set_outputMatrix(obj,C)
            obj.C = C;
        end
        
        % set uTrans
        function set_uTrans(obj, uTrans, loc)
            obj.options.uLocTrans{loc} = uTrans;
        end
        
        % set input set U for specific location
        function set_Uloc(obj, U, loc)
            obj.options.Uloc{loc} = U;
        end
        
        % set input set U globally (same for each location)
        function set_U(obj, U)
            for i=1:obj.modes
                obj.set_Uloc(U,i)
            end
        end
        
        % set guardIntersect
        function set_guarIntersect(obj,guardIntersect)
            obj.options.guardIntersect = guardIntersect;
        end
        
        % set isHyperplaneMap
        function set_HyperplaneMap(obj,num)
            obj.isHyperplaneMap = num;
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
        
        % set originContained
        function set_originContained(obj, originContained)
            obj.options.originContained = originContained;
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
        
        % set inital state for simulation x0
        function set_x0(obj, x0)
            if length(x0) ~= obj.dim
                error('Dimension mismatch between initial state and the system');
            end
            obj.options.x0 = x0;
        end
        
        % set time step 
        function set_timeStep(obj, timeStep)
            if timeStep <= 0
                error('Invalid time step');
            end
            obj.options.timeStep = timeStep;
        end
        
    end
    
    
    % reachability anlaysis method
    methods
        
        % reachability analysis using zonotope
        function [R, reachTime] = reach_zono(obj, init_set, input_set, timeStep, tFinal)
           % [R, reachTime] = reach_zono(obj, init_set, input_set, timeStep, tFinal)
           % @init_set: initial set of state
           % @input_set: input set u
           % @timeStep: time step in reachable set computation
           % @num_steps: number of steps in reachable set computation
           % R: output reachable set
           % reachTime = time duration of reachability analysis
           %
           % Example
           % [R, rT] = HA.reach_zono(init_set, input_set, timeStep, tFinal);
           
           start = tic;
           if ~isa(init_set, 'Zono') 
               error('Initial set is not a zonotope');
           end
           
           if ~isa(input_set, 'Zono')
               error('Input set is not a zonotope');
           end
           
           R0 = zonotope([init_set.c, init_set.V]);
           if input_set.dim == 0
               U = zonotope(0);
           else
            U = zonotope([input_set.c, input_set.V]);
           end
           
           obj.set_R0(R0);
           obj.set_U(U);
           obj.set_timeStep(timeStep);
           obj.set_tFinal(tFinal);
           
           R = reach(obj.sysCORA, obj.options); % CORA reach method using zonotope and conservative linearization
           R = get(R,'reachableSet');
           R = R.OT;
           reachTime = toc(start);
        end
        
        % step reach using star set used for neural network control system
        function S = stepReachStar(obj, init_set, input_set)
            % S = stepReachStar(obj, init_set, input_set)
            % @init_set: initial set, a star
            % @input_set: input set, a star
            % S: reachable set at the timeStep, a star
            
            if ~isa(init_set, 'Star')
                error('Initial set is not a star');
            end
            
            if ~isa(input_set, 'Star')
                error('Input set is not a star');
            end
            
            I = init_set.getZono;
            if input_set.dim == 0
                U = Zono;
            else
                U = input_set.getZono;
            end
            
            [R, ~] = obj.reach_zono(I, U, obj.options.timeStep, obj.options.tFinal);
            
            N = length(R);
            Z = R{N}{1,1}; % get the last zonotope in the reach set
            Z = Z.Z; % get c and V 
            c = Z(:,1); % center vector
            V = Z(:, 2:size(Z, 2)); % generators
            
            Z = Zono(c, V);
            S = Z.toStar;
            
            for i=1:N
                N = length(R);
                Z = R{i}{1,1}; 
                Z = Z.Z; % get c and V 
                c = Z(:,1); % center vector
                V = Z(:, 2:size(Z, 2)); % generators

                Z = Zono(c, V);
                S = Z.toStar;
                obj.intermediate_reachSet = [obj.intermediate_reachSet S];
            end
            % the last zonotope in the reach set is returned
        end
        
        % evaluate (simulate) the plant
        function y = evaluate(obj,u, x0)
            % y = evaluate(obj, u, x0)
            % simulate hybrid automaton
            obj.set_x0(x0);
            for i=1:obj.modes
                obj.options.uLoc{i} = u;
            end
            HAsim = simulate(obj.sysCORA,obj.options);
            y = get(HAsim,'trajectory');
        end
        
    end
end

