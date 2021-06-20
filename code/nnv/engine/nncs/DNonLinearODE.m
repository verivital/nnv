classdef DNonLinearODE < handle
    % Discrete time nonlinear ODE class. Grapper of nonlinearSysDT in CORA.
    % Chelsea Sidrane April 25, 2019 && Dung Tran: 4/29/2019
    % Last Revision: September 3, 2020 - Diego Manzanas
    %      - CORA 2020 updates
    % Based on NonLinearODE.m
    
    properties
        options = []; % option for recahable set computation
        params = []; 
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
                        
            % default options and parameters
            param.tStart = 0;
            param.tFinal = Ts; % we provide method to change the option
            param.R0 = []; % initial set for reachability analysis
            option.zonotopeOrder = 2; % zonotope order
            option.reductionTechnique = 'girard';
            option.errorOrder = 1;
            option.tensorOrder=2;
            % Store options and parameters
            obj.options = option; % default option
            obj.params = param; % default parameters
            obj.set_output_mat(outputMat);  
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
        function [R, reachTime] = reach_zono(obj, init_set, input_set)
            % Syntax:
            % [R, reachTime] = reach_zono(obj, init_set, input_set)
            %
            % Inputs:
           % init_set: initial set of state
           % input_set: input set u
           %
           % Outputs: 
           % R = reachable set of the states of the plant
           % reachTime = total time during reachability computation
           
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
           dt  = obj.Ts;
           fun = @(x,u) obj.dynamics_func(x,u,dt);
           sys = nonlinearSysDT(fun, obj.Ts, obj.dim, obj.nI); % CORA nonlinearSys class
           R = reach(sys, obj.params, obj.options); % CORA reach method using zonotope and conservative linearization
                     
           reachTime = toc(start);
                   
        end
        
       % step reach using star set used for neural network control system
        function S = stepReachStar(obj, init_set, input_set)
            % Syntax:
            % S = stepReachStar(obj, init_set, input_set)
            % 
            % Inputs:
            % @init_set: initial set, a star
            % @input_set: input set, a star
            %
            % Outputs:
            % @R: reachable set, a star
            
            if ~isa(init_set, 'Star')
                error('Initial set is not a star');
            end
            
            if ~isa(input_set, 'Star')
                error('Input set is not a star');
            end
            
            I = init_set.getZono;
            U = input_set.getZono;
            
            [R, ~] = obj.reach_zono(I, U);
            
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
                Z = R(i).timePoint.set; % get the reachset 
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
        end
        % evaluate (simulate) the plant with specific input and state
        % using ode45 solver
        function y = evaluate(obj, x0, u)
            % Syntax:
            % y = evaluate(obj, x0, u)
            %
            % Inputs:
            % @x0: initial state
            % @u: control input
            % 
            % Outputs:
            % y: state trajectory
            
            % author: Dung Tran
            % date: 1/29/2019
            
            % y = obj.dynamics_func(x0, u);
            y = obj.dynamics_func(x0, u);
            
        end
        
        % implement box?
        
        % implement polyhedron?
        
        % add simulation based evaluation?
        
    end
end