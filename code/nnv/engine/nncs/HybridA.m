classdef HybridA < handle
    %HYBRID AUTOMATA class is a wrapper of Hybrid Automaton class in CORA 
    %   Author: Diego Manzanas 05/07/2020
    %   Last revision: Diego Manzanas 09/11/2020
    % 
    %HYBRIDA Construct an instance of this class
    %      HA = HybridA(dum, nI,dynamics_func, modes);
    %      HA = NNV HA wrapper class for CORA Hybrid Automaton class
    %      dim = dynamics dimensions for each mode
    %      nI = number of inputs
    %      dynamics_func = hybrid dynamics in CORA format
    %      modes = number of modes of the hybrid automata dynamics
    
    properties
        options = []; % option for reachable set computation
        params = []; % paramters for reachable set computations
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
            %HYBRIDA Construct an instance of this class
            %   HA = HybridA(dim, nI, dynamics_func, modes);
            %   HA = NNV HA wrapper class for CORA Hybrid Automaton class
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
        
        % Initialize the HybridA reachability options to default (some are
        % not used in different examples, could initialize them depending
        % if linear or nonlinear system)
        % a = get(HA.sysCORA, 'location');
        % b = get(a,'contDynamics');
        % if contains(class(b),'non') -> nonlinear sys options
        function HybridA_init(obj)
            % Reachability Time default option
            param.tStart = 0;
            param.tFinal = 0.1; % we provide method to change the option
            option.timeStep = 0.05;
            % Reachability mode location default option
            param.startLoc = 1; %initial location
            param.finalLoc = 0; %0: no final location
            
%             param.x0 = []; %initial state for simulation
            param.R0 = []; %initial state for reachability analysis
            for i=1:obj.modes
%                 option.timeStepLoc{i} = 0.05; %time step size for reachable set computation in location i
                % Initialize inputs (no inputs)
%                 param.uLoc{i} = 0;
                param.uLocTrans{i} = 0;
%                 param.Uloc{i} = zonotope(0);
            end
            option.taylorTerms = 10;
            option.zonotopeOrder = 20;
            option.errorOrder = 2;
            option.reductionTechnique = 'girard';
            option.guardIntersect = 'polytope';
            option.linAlg = 'standard';
            option.intermediateOrder = 2;
            option.tensorOrder = 2;
            option.reductionInterval = inf;
            option.maxError = 1e7*ones(obj.dim,1);
            option.enclose = {'box'};
            option.alg = 'lin';
            option.guardOrder = 2;
            option.intersectInvariant = 0;
            obj.options = option; % default options
            obj.params = param; % default parameters
        end
        
        % set advanced linear error comp.
        function set_advancedLinErrorComp(obj,aLE)
            obj.options.advancedLinErrorComp = aLE;
        end
        
        % set reduction interval
        function set_reductionInterval(obj,rI)
            obj.options.reductionInterval = rI;
        end
        
        % set old error
        function set_oldError(obj,oE)
           obj.options.oldError = oE; 
        end
        
        % set maximum error
        function set_maxError(obj,mE)
           obj.options.maxError = mE; 
        end
        
        % set tensor order
        function set_tensorOrder(obj,tO)
           obj.options.tensorOrder = tO; 
        end
        
        % set intermediate order
        function set_intermediateOrder(obj,iO)
           obj.options.intermediateOrder = iO; 
        end
        
        % set linAlg
        function set_linAlg(obj,a)
            obj.options.linAlg = a;
        end
        
        % set max projection error
        function set_maxProjectionError(obj,a)
            obj.maxProjectionError = a;
        end

        % set enclonsure enable method
        function set_enclosureMethod(obj,num)
            obj.options.enclosureEnables = num;
        end
        
        % set output matrix
        function set_outputMatrix(obj,C)
            obj.C = C;
        end
        
        % set uTrans
        function set_uTrans(obj, uTrans, loc)
            obj.params.uLocTrans{loc} = uTrans;
        end
        
        % set input set U globally (same for each location)
        function set_U(obj, U)
            obj.params.U = U;
        end
        
        % set guardIntersect
        function set_guarIntersect(obj,guardIntersect)
            obj.options.guardIntersect = guardIntersect;
        end
        
        % set isHyperplaneMap
        function set_HyperplaneMap(obj,num)
            obj.options.isHyperplaneMap = num;
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
        
        % set start location
        function set_startLoc(obj,loc)
            obj.params.startLoc = loc;
        end
        
        % set final location
        function set_finalLoc(obj,loc)
            obj.params.finalLoc = loc;
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
           % @tFinal: final time of reachability computation
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
           
           R = reach(obj.sysCORA, obj.params, obj.options); % CORA reach method using zonotope and conservative linearization
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
            
            [R, ~] = obj.reach_zono(I, U, obj.options.timeStep, obj.params.tFinal);
            N = length(R); % number of reachsets in the computation
            max_parent = 0;
            for pp=1:N
                max_parent = max(max_parent,R(pp).parent);
            end
            S = [];
            for mp=1:N
                if R(mp).parent == max_parent
                    Z = R(mp).timePoint.set; % get the last reachset
                    Nn = length(Z); % number of sets in the last reachset
                    Z = Z{Nn}; % get the last set in the reachset
                    Nz = length(Z); % number of zonotopes in the last set
                    for ik=1:Nz
                        try
                            Z1 = Z{ik}.Z;
                        catch
                            Z1 = Z.Z; % get c and V 
                        end
                        c = Z1(:,1); % center vector
                        V = Z1(:, 2:size(Z1, 2)); % generators

                        Zz = Zono(c, V);
                        S = [S Zz.toStar];
                    end
                end
            end
            
            for i=1:N
%                 N = length(R);
                Z = R(i).timeInterval.set; % get the reachset 
%                 Z = R(i).timePoint.set;
                Nn = length(Z); % number of sets in the reachset (1 x timeStep)
                Ss = [];
                for ik=1:Nn
                    Z1 = Z{ik};
                    Nz = length(Z1);
                    for iz=1:Nz
                        try
                            Z2 = Z1{iz}.Z;
                        catch
                            Z2 = Z1.Z; % get c and V 
                        end
    %                     Z = Z.Z; % get c and V 
                        c = Z2(:,1); % center vector
                        V = Z2(:, 2:size(Z2, 2)); % generators

                        Zz = Zono(c, V);
                        Ss = [Ss Zz.toStar];
                    end
                end
                obj.intermediate_reachSet = [obj.intermediate_reachSet Ss];
            end
            
            % the last zonotope in the reach set is returned
            
        end
        
        % evaluate (simulate) the plant
        function [t,x,loc] = evaluate(obj,u, x0)
            % [t,x,loc] = evaluate(obj, u, x0)
            % simulate hybrid automaton
            simOpt = obj.options;
            simOpt.x0 = x0;
            for i=1:obj.modes
                simOpt.uLocTrans{i} = u;
            end
            simOpt.tStart = obj.params.tStart;
            simOpt.tFinal = obj.params.tFinal;
            simOpt.startLoc = obj.params.startLoc;
            simOpt.u = u;
            [t,x,loc] = simulate(obj.sysCORA, simOpt);
        end
        
    end
end

