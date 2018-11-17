classdef NNCS < handle
    %Neural network control system class 
    %   Dung Tran: 10/21/2018
    
    properties
        controller = []; % nerual-network controller
        plant = []; % plant model, could be linear, nonlinear or neural network-based plant
        feedbackMap = []; % a feedback matrix decribes the mapping from a group of 
                          % outputs of the plant to a group of inputs of the controller
                          
        % nerual network control system architecture
        %
        %              ---> plant ---> y(t) ---sampling--->y(k) 
        %             |                                       |
        %             |                                       |
        %             u(k) <---- controller |<---- y(k-d)-----(output feedback) 
        %                                   |<----- v(k)------(reference input)                                    
        
        
        % the input to neural net controller is grouped into 2 group
        % the first group contains all the reference inputs
        % the second group contains all the output feedback with delays
        
        % feedbackMap = [0;1], a 2 x 1 matrix, means that:
        % the output feedback to the controller are: [y[k]; y[k-1]]        
        
        % the first layer weight matrix of the controller is decomposed into two
        % submatrices: W = [W1 W2] where
        %              W1 is conresponding to I1 = v[k] (the reference input)
        %              W2 is conresponding to I2 = [y[k]; y[k-1]] (the feedback inputs)  
        
        % the reach set of the first layer of the controller is: 
        %              R = f(W1 * I1 + W2 * I2 + b), b is the bias vector of
        %              the first layer, f is the activation function
        
        nO = 0; % number of output
        nI = 0; % number of inputs = size(I1, 1) + size(I2, 1)
        nI_ref = 0; % number of reference inputs
        nI_fb = 0; % number of feedback inputs
        
        % used for reachable set computation
        ref_I = []; % reference input set
        init_set = []; % initial set for the plant
        reachSetTree = []; % reachable set tree
        totalNumOfReachSet = 0; % total number of reachable sets
        reachTime = 0; % reachable set computation time
    end
    
    methods
        
        %constructor
        function obj = NNCS(controller, plant, feedbackMap)
            % @controller: a neural net controller
            % @plant: a plant model (a LinearODE, DLinearODE or Neural net)
            % @feedbackMap: a feedback map from outputs of the plant to the
            % input of the controller
            
            % author: Dung Tran
            % date: 11/1/2018
            
            
            if ~isa(controller, 'FFNN')
                error('The controller is not a feedforward neural network');
            end
            
            if ~isa(plant, 'LinearODE') && ~isa(plant, 'DLinearODE')
                error('The plant is not a linear ode system or a discrete linear ode system');
            end            
                        
            [nO, nI] = size(feedbackMap);
            
            if nI ~= 1
                error('FeedbackMap should have one column');
            end
            if nO * plant.nO > controller.nI
                error('Two many feedback inputs');
            end
                        
            obj.controller = controller;
            obj.plant = plant;
            obj.feedbackMap = feedbackMap;
            obj.nO = plant.nO;
            obj.nI = controller.nI;
            obj.nI_fb = nO * plant.nO;
            obj.nI_ref = controller.nI - obj.nI_fb;
            
        end
        
         % reachability analysis of NNCS using polyhedron
         % main limitation: the number of vertices in the reachable set
         % increase quickly which affects the scalability of this approach
        function [P, reachTime] = reachPolyhedron_approx(obj, init_set, ref_inputSet, n_steps)
             % @init_set: the initial set of condition for the plant
             % @ref_inputSet: the reference input set applied to the controller
             % @n_steps: number of steps 
             % @P: the state reachable set of the plant
             %     we get the output reachable set by mapping P on the
             %     direction of interest plant.C

             % author: Dung Tran
             % date: 11/2/2018

             start_time = tic; 
             if ~isa(obj.plant, 'DLinearODE')
                 error('Reachability analysis of NNCS using Polyhedron only supports for Discrete linear ODE plant');
             end

             if ~isa(init_set, 'Polyhedron')
                 error('Initial set of the plant is not a polyhedron');
             end

             if ~isempty(ref_inputSet) && ~isa(ref_inputSet, 'Polyhedron')
                 error('The reference input set is not a polyhedron');
             end

             if n_steps < 1
                 error('Number of steps should be >= 1');
             end
             
             obj.reachSetTree = SetTree(n_steps + 1); % initialize reach set tree
             obj.init_set = init_set;
             obj.ref_I = ref_inputSet;
             obj.reachSetTree.addReachSet(init_set, 1); % add the init_set to the reachSetTree

             for i=2:n_steps + 1                 
                 fb_I = obj.reachSetTree.extract_fb_ReachSet(i - 1);   
                 input_set = obj.nextInputSetPolyhedron(fb_I{1});
                 [U,~] = obj.controller.reach(input_set, 'exact', 1, []); % control set at step i
                 U1 = Reduction.fastHull(U);
                 %U1 = Reduction.hypercubeHull(U);
                 %U1 = U1.toPolyhedron;
                 R = obj.plant.stepReachPolyhedron(fb_I{1}(length(fb_I{1})), U1);                 
                 obj.reachSetTree.addReachSet(R, i);                 
             end
             reachTime = toc(start_time);
             obj.reachTime = reachTime;
             P = obj.reachSetTree.flatten();
             obj.totalNumOfReachSet = obj.reachSetTree.getTotalNumOfReachSet();
             
        end
        
        
        % reachability analysis of nncs using zonotope
        % output reach set of controller is a single box
        % the plant reachable set is a zonotope
        function [P, reachTime] = reach_zono(obj, init_set, ref_inputSet, n_cores, n_steps)
             % @init_set: the initial set of condition for the plant
             % @ref_inputSet: the reference input set applied to the controller
             % @n_steps: number of steps 
             % @P: the state reachable set of the plant
             %     we get the output reachable set by mapping P on the
             %     direction of interest plant.C

             % author: Dung Tran
             % date: 11/16/2018
             
             start_time = tic; 
             if ~isa(obj.plant, 'DLinearODE')
                 error('Reachability analysis of NNCS using Polyhedron only supports for Discrete linear ODE plant');
             end

             if ~isa(init_set, 'Polyhedron')
                 error('Initial set of the plant is not a Box');
             end

             if ~isempty(ref_inputSet) && ~isa(ref_inputSet, 'Polyhedron')
                 error('The reference input set is not a polyhedron');
             end

             if n_steps < 1
                 error('Number of steps should be >= 1');
             end
             
             obj.reachSetTree = SetTree(n_steps + 1); % initialize reach set tree
             obj.init_set = init_set;
             obj.ref_I = ref_inputSet;
             obj.reachSetTree.addReachSet(init_set, 1); % add the init_set to the reachSetTree

             for i=2:n_steps + 1                 
                 fb_I = obj.reachSetTree.extract_fb_ReachSet(i - 1);   
                 input_set = obj.nextInputSetPolyhedron(fb_I{1});
                 [U,~] = obj.controller.reach(input_set, 'exact', 1, []); % control set at step i
                 U1 = Reduction.hypercubeHull(U);
                 % to do:
                 plant_I = fb_I{1}(length(fb_I{1}));
                 R = obj.plant.stepReachBox(plant_I, U1);                 
                 obj.reachSetTree.addReachSet(R, i);                 
             end
             reachTime = toc(start_time);
             obj.reachTime = reachTime;
             P = obj.reachSetTree.flatten();
             obj.totalNumOfReachSet = obj.reachSetTree.getTotalNumOfReachSet();
             
             
        end
        
        
        % reach Polyhedron with exact reachable set computation
        
        function [P, reachTime] = reachPolyhedron_exact(obj, init_set, ref_inputSet, n_cores, n_steps)
             % @init_set: the initial set of condition for the plant
             % @ref_inputSet: the reference input set applied to the controller
             % @n_cores: number of cores used in computation
             % @n_steps: number of steps
             % @Px: the state reachable set of the plant
             % @Py: the output reachable set of the plant
             %     we get the output reachable set by mapping Px on the
             %     dimension of interest

             % author: Dung Tran
             % date: 11/6/2018
            
             start_time = tic; 
             
             if ~isa(obj.plant, 'DLinearODE')
                 error('Reachability analysis of NNCS using Polyhedron only supports for Discrete linear ODE plant');
             end

             if ~isa(init_set, 'Polyhedron')
                 error('Initial set of the plant is not a polyhedron');
             end

             if ~isempty(ref_inputSet) && ~isa(ref_inputSet, 'Polyhedron')
                 error('The reference input set is not a polyhedron');
             end
             
             if n_steps < 1
                 error('number of step should be >= 1');
             end
             
             if n_cores < 1
                 error('number of cores should be >= 1');
             end
             
             obj.reachSetTree = SetTree(n_steps + 1); % initialize reach set tree
             obj.init_set = init_set;
             obj.ref_I = ref_inputSet;
             
             obj.reachSetTree.addReachSet(init_set, 1); % add the init_set to the reachSetTree
             
             for i=2:n_steps + 1
                 obj.stepReachPolyhedron(i, n_cores);                 
             end
             
             reachTime = toc(start_time);
             obj.reachTime = reachTime;
             P = obj.reachSetTree.flatten();
             obj.totalNumOfReachSet = obj.reachSetTree.getTotalNumOfReachSet();
             
        end
        
        
        % step reach Polyhedron
        function R = stepReachPolyhedron(obj, i, n_cores)
            % @i: the step number
            % @n_cores: number of cores used in computation
            % @R: the reachable set of step i   
            
            % author: Dung Tran
            % date: 11/6/2018
            
            if i < 1
                error('Step ID should >= 1');
            end
            
            fb_R = obj.reachSetTree.extract_fb_ReachSet(i - 1); % feedback reach set 
            % the last element of fb_R{1, i} is the initial set of the
            % plant for the current step
            n = length(fb_R);
            
            R = [];            
            for j=1:n
                l = length(fb_R{1, j});
                cur_init_set = fb_R{1, j}(l); % init set for current step
                R1 = obj.stepReachPolyhedron_singleFeedBack(cur_init_set, fb_R{1, j}, n_cores);
                R = [R R1];               
            end
                        
            obj.reachSetTree.addReachSet(R, i);        
            
        end
     
        % step reach Polyhedron with a single feedback set
        function R = stepReachPolyhedron_singleFeedBack(obj, init_set, fb_I, n_cores)
            % @init_set: the initial set of the step
            % @ref_I: the reference input set
            % @fb_I_cell: the feedback input set
            % @n_cores: number of cores used for computation
            % @R:  state reachable set
            % @fb_R: the feedback input set for the next step
            
            % author: Dung Tran
            % date: 11/6/2018
            
            
            input_set = obj.nextInputSetPolyhedron(fb_I); % get the input set for the current step
            
            if isempty(input_set) && isempty(init_set)
                R = [];
            end
            
            if isempty(input_set) && ~isempty(init_set)
                R = obj.plant.stepReachPolyhedron(init_set, []);             
            end
            
            if ~isempty(input_set)
                % compute reachable set for the controller
                U = obj.controller.reach(input_set, 'exact', n_cores, []);
                % U is an array of polyhedra which are the input sets to the
                % plant
                
                % compute reachable set for the plant
                R = obj.plant.stepReachPolyhedron_parallel(init_set, U, n_cores);
                
            end
      

            
        end
        
               
        % get the next step input set
        function I = nextInputSetPolyhedron(obj, fb_I)
            % @fb_I: the feedback input set
            
            % author: Dung Tran
            % date: 11/5/2018
            
            l = length(fb_I);
            fb_inputSet = [];
            if l > 0
                for i=1:l
                    fb_inputSet = [fb_inputSet fb_I(i).affineMap(obj.plant.C, 'vrep')];
                end
            end
            n = size(obj.feedbackMap, 1);          
            
            if isempty(obj.ref_I) && isempty(fb_inputSet)
                I = [];
            end
            if ~isempty(obj.ref_I) && isempty(fb_inputSet)
                
                lb = zeros(obj.nI_fb);
                ub = zeros(obj.nI_fb);
                
                I2 = Polyhedron('lb', lb, 'ub', ub);
                I = Conversion.concatenatePolyhedron([obj.ref_I I2]);
                
            end
            
            if ~isempty(obj.ref_I) && ~isempty(fb_inputSet)
               
                l = length(fb_inputSet);
                nA = size(fb_inputSet(1).A,2);
                I2 = [];
                for i=1:n

                    if l - obj.feedbackMap(i) <= 0
                        I2 = [I2 Polyhedron('lb', zeros(nA, 1), 'ub', zeros(nA, 1))];
                    else

                        I2 = [I2 fb_inputSet(l - obj.feedbackMap(i))];

                    end                

                end

                I = Conversion.concatenatePolyhedron([obj.ref_I I2]);
            end
            
            
            if isempty(obj.ref_I) && ~isempty(fb_inputSet)
                
                l = length(fb_inputSet);
                nA = size(fb_inputSet(1).A,2);
                I2 = [];
                for i=1:n

                    if l - obj.feedbackMap(i) <= 0
                        I2 = [I2 Polyhedron('lb', zeros(nA, 1), 'ub', zeros(nA, 1))];
                    else

                        I2 = [I2 fb_inputSet(l - obj.feedbackMap(i))];

                    end                

                end
                
                lb = zeros(obj.nI_ref,1);
                ub = zeros(obj.nI_ref,1);
                I1 = Polyhedron('lb', lb, 'ub', ub);
                I = Conversion.concatenatePolyhedron([I1 I2]);
                
            end
            
        end
                              
        
    end
    
    
    methods(Static)
        
    end
    
end

