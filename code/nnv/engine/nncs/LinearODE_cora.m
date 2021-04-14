classdef LinearODE_cora < handle
    % Linear system represented by ODE x' = Ax + Bu, y = Cx + Du
    %   Linear ODE (CORA methods) system class
    %
    %
    % Construct LinearODE plant (CORA)
    % 
    % Syntax:
    %     obj = LinearODE_cora(A,B,C,D,reachTimeStep, controlPeriod)
    %
    %     @A: system matrix A
    %     @B: system matrix B
    %     @C: dynamics of the plant, the input should have
    %     @D: character to specify function, for example @car_dynamics
    %     @reachTimeStep: reachability step for the plant
    %     @controlPeriod: control period for the plant
    %
    % Example:
    %     A = [-2 0; 1 -3];
    %     B = [1; 1];
    %     C = [1 0];
    %     D = [0];
    %     LinearPlant = LinearODE_cora(A,B,C,D,0.05,0.2);
    %
    % author: Diego Manzanas
    % date: 10/12/2020

    % Note: we construct a LinearODE plant with default option (CORA)
    
    properties
        A = []; % system matrix A
        B = []; % system matrix B
        C = []; % system matrix C
        D = []; % system matrix D
        dim = 0; % system dimension
        nI = 0; % number of inputs
        nO = 0; % number of outputs
        controlPeriod = 0.1; % control period
        reachStep = 0.01; % reachability step for the plant
        intermediate_reachSet = []; % Store all computed reach sets
        options = []; % reachability options (CORA)
        params = []; % reachability parameters (CORA)
    end
    
    methods
         % constructor from a matlab function reference in dynamics_func
        function obj = LinearODE_cora(A,B,C,D, reachTimeStep, controlPeriod)
            % construct LinearODE plant (CORA)
            % 
            % Syntax:
            % obj = LinearODE_cora(A,B,C,D,reachTimeStep, controlPeriod)
            %
            % @A: system matrix A
            % @B: system matrix B
            % @C: system matrix C
            % @D: system matrix D
            % @reachTimeStep: reachability step for the plant
            % @controlPeriod: control period for the plant
            %
            % Example:
            % A = [-2 0; 1 -3];
            % B = [1; 1];
            % C = [1 0];
            % D = [0];
            % LinearPlant = NonLinearODE_cora(A,B,C,D,0.05,0.2);
            %
            % author: Diego Manzanas
            % date: 10/12/2020
            
            % Note: we construct a LinearODE plant with default option (CORA)
            
            if isempty(A)
                error('Matrix A is empty');
            else 
                [nA, mA] = size(A);
                if nA ~= mA
                    error('A should be a square matrix');
                end
            end
            
            if ~isempty(B)
                [nB, ~] = size(B);
                if nA ~= nB
                    error('A and B are inconsistent');
                end
            end
            
            if ~isempty(C)
                [nC, mC] = size(C);
                if mC ~= nA
                    error('A and C are inconsistent');
                end
            end
            
            if ~isempty(C) && ~isempty(D)
               [nD, ~] = size(D); 

                if nD ~= nC
                    error('C and D are inconsistent');
                end
            end
            
            obj.A = A;
            obj.B = B;
            obj.C = C;
            obj.D = D;
            
            obj.LinearODE_init();
            obj.set_timeStep(reachTimeStep);
            obj.reachStep = reachTimeStep;
            obj.set_tFinal(controlPeriod);
            obj.controlPeriod = controlPeriod;
        end
        
        % set default parameters
        function LinearODE_init(obj)
            % default option
            param.tStart = 0;
            param.tFinal = 0.2; % we provide method to change the option
            
            param.R0 = []; % initial set for reachability analysis
            param.U = []; % input set for reachability analysis
            option.timeStep = 0.01; % time step for reachable set computation
            option.taylorTerms = 4; % number of taylor terms for reachable sets
            option.zonotopeOrder = 10; % zonotope order
            option.reductionTechnique = 'girard';
            option.linAlg = 'standard';
            
            obj.options = option; % default option
            obj.params = param; % default parameters
        end
        
        % Set parameters and options
        function set_timeStep(obj, tStep)
            obj.options.timeStep = tStep;
        end
        
        function set_tFinal(obj,tFinal)
            obj.params.tFinal = tFinal;
        end
        
        function set_zonotopeOrder(obj,order)
            obj.options.zonotopeOrder = order;
        end
        
        function set_taylorTerms(obj,terms)
            obj.options.taylorTerms = terms;
        end
        
        function set_reductionTechnique(obj,technique)
            obj.options.reductionTechnique = technique;
        end
        
        function set_linAlg(obj,alg)
            obj.options.linAlg = alg;
        end
        
        function set_R0(obj,R)
            obj.params.R0 = R;
        end
        
        function set_U(obj,U)
            obj.params.U = U;
        end
    end
    
    % reachability anlaysis method
    methods
        
        % reachability analysis using zonotope
        function [R, reachTime] = reach_zono(obj, init_set, input_set, timeStep, tFinal)
           % @init_set: initial set of state
           % @input_set: input set u
           % @timeStep: time step in reachable set computation
           % @tFinal:final time in reachable set computation
           
           % this is a grapper of reach method for linear system in CORA
           % the initial set and input set are Zono in nnv
           
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
           
           sys = linearSys(obj.A, obj.B, [], obj.C, obj.D); % CORA linearSys class
           R = reach(sys, obj.params, obj.options); % CORA reach method using zonotope and conservative linearization
                     
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
%                 Z = R(i).timeInterval.set; % get the reachset 
                Z = R(i).timePoint.set;
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
            sys = linearSys(obj.A, obj.B, [], obj.C, obj.D);
            [t,y] = simulate(sys, simOpt);
        end
    end
end
