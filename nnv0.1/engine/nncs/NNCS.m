classdef NNCS
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
        function [Px, Py] = reachPolyhedron(obj, init_set, ref_inputSet, N)
             % @init_set: the initial set of condition for the plant
             % @ref_inputSet: the reference input set applied to the controller
             % @N: number of steps 
             % @Px: the state reachable set of the plant
             % @Py: the output reachable set of the plant
             %     we get the output reachable set by mapping Px on the
             %     dimension of interest

             % author: Dung Tran
             % date: 11/2/2018


             if ~isa(obj.plant, 'DLinearODE')
                 error('Reachability analysis of NNCS using Polyhedron only supports for Discrete linear ODE plant');
             end

             if ~isa(init_set, 'Polyhedron')
                 error('Initial set of the plant is not a polyhedron');
             end

             if ~isempty(ref_inputSet) && ~isa(ref_inputSet, 'Polyhedron')
                 error('The reference input set is not a polyhedron');
             end

             if N == 0
                 P = init_set;
             end

             if N == 1
                 Rx = obj.plant.stepReachPolyhedron(init_set, []); % state reachable set
                 Px = [init_set Rx];
                 Py = [init_set.affineMap(obj.plant.C) Rx.affineMap(obj.plant.C)];
             end         

             if N > 1

                 % first step
                 Rx = obj.plant.stepReachPolyhedron(init_set, []); % state reachable set
                 Px = [init_set Rx];
                 Py = [init_set.affineMap(obj.plant.C) Rx.affineMap(obj.plant.C)];


                 for i=2:N
                     % next input set for the controller
                     I = obj.nextInputSetPolyhedron(ref_inputSet, Py); % input set to the controller at step i
                     U = obj.controller.reach(I, 'exact', 1, []); % control set at step i
                     U1 = Reduction.fastHull(U);
                     Rx = obj.plant.stepReachPolyhedron(Px(length(Px)), U1); 
                     Px = [Px Rx];
                     Ry = Rx.affineMap(obj.plant.C);
                     Py = [Py Ry];                 
                 end

             end


        end
        
        % reach Polyhedron with exact reachable set computation
        
        function [Px, Py] = reachPolyhedron_exact(obj, init_set, ref_inputSet, n_cores, n_steps)
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
            
             if ~isa(obj.plant, 'DLinearODE')
                 error('Reachability analysis of NNCS using Polyhedron only supports for Discrete linear ODE plant');
             end

             if ~isa(init_set, 'Polyhedron')
                 error('Initial set of the plant is not a polyhedron');
             end

             if ~isempty(ref_inputSet) && ~isa(ref_inputSet, 'Polyhedron')
                 error('The reference input set is not a polyhedron');
             end

             if n_steps == 0
                 Px = init_set;
                 Py = init_set.affineMap(obj.plant.C, 'vrep');
             end

             if n_steps > 0
                 
                 % step 0
                 Px = init_set;
                 
                 for i=1:n_steps
                     
                    if i == 1
                        % step 1: reference_input = []; feedback_input = []
                        ref_I = [];
                        fb_I = cell(1,1);
                        [R, fb_R] = obj.stepReachPolyhedron_singleFeedBack(init_set, ref_I, fb_I_cell, n_cores);
                        Px = [Px R];                     
                    else
                        
                        n = length(R); % multiple init_set for the plant
                        m = length(fb_R); % m feedback set for the controller
                        
                        % n init_set combined with m feedback set for the
                        % controller produce p reachable set, p is unknown
                        R2 = [];
                        fb2 = [];
                        for j=1:n                            
                            [R1, fb1] = obj.stepReachPolyhedron_multipleFeedbacks(R(j), ref_I, fb_R, n_cores);                            
                            R2 = [R2 R1];
                            fb2 = [fb2 fb1];                            
                        end
                        
                        R = R2;
                        fb_R = fb2;
                        Px = [Px R];
                        
                    end

                     
                 end
                 
                 
                 
             end
             


             
             
             
            
        end
     
        % step reach Polyhedron with a single feedback set
        function [R, fb_R] = stepReachPolyhedron_singleFeedBack(obj, init_set, ref_I, fb_I_cell, n_cores)
            % @init_set: the initial set of the step
            % @ref_I: the reference input set
            % @fb_I_cell: the feedback input set
            % @n_cores: number of cores used for computation
            % @R:  state reachable set
            % @fb_R: the feedback input set for the next step
            
            % author: Dung Tran
            % date: 11/6/2018
            
            fb_I = fb_I_cell{1, 1}; % get an array of feedback set 
            
            input_set = obj.nextInputSetPolyhedron(ref_I, fb_I); % get the input set for the current step
            
            if isempty(input_set) && isempty(init_set)
                R = [];
                fb_R = [];
            end
            
            if isempty(input_set) && ~isempty(init_set)
                R = obj.plant.stepReachPolyhedron(init_set, []);
                fb_R = R.affineMap(obj.plant.C, 'vrep');             
            end
            
            if ~isempty(input_set)
                % compute reachable set for the controller
                U = obj.controller.reach(input_set, 'exact', n_cores, []);
                % U is an array of polyhedra which are the input sets to the
                % plant
                
                % compute reachable set for the plant
                p = length(U);
                fb_R = cell(1, p);
                R = obj.plant.stepReachPolyhedron_parallel(init_set, U, n_cores);
                for i=1:length(R)
                    fb_R{1, i} = [fb_I R(i).affineMap(obj.plant.C, 'vrep')];
                end
                
            end
      

            
        end
        
        % step Reach Polyhedron with multiple feedback sets
        function [R, fb_R] = stepReachPolyhedron_multipleFeedbacks(obj, init_set, ref_I, fb_Is, n_cores)
            % @init_set: the initial set of plant
            % @ref_I: the reference input
            % @fb_Is: a cell of feedback sets
            % @n_cores: number of cores used in computation
            % @R: the state reachable set
            % @fb_R: the output feedback set for the next step
            
            % author: Dung Tran
            % date: 11/8/2018
            
            n = length(fb_Is); % an array of cell
            if n == 0
                [R, fb_R] = obj.stepReachPolyhedron_singleFeedback(init_set, ref_I, [], n_cores);
            else
                fb_R = []; % an array of cell
                R = [];
                for i=1:n                  
                   [R1, fb_R1] = obj.stepReachPolyhedron_singleFeedback(init_set, ref_I, fb_Is(i), n_cores);
                   R = [R R1];
                   fb_R = [fb_R fb_R1];
                end
                
            end
            
            
            
        end
     
        % get the next step input set
        function I = nextInputSetPolyhedron(obj, ref_I, outputSet)
            % @ref_I: is the reference input
            % @outputSet: the output set
            % initial step ref_I = [] and outputSet = []
            
            % author: Dung Tran
            % date: 11/5/2018
            
            n = size(obj.feedbackMap, 1);          
            
            if isempty(ref_I) && isempty(outputSet)
                I = [];
            end
            if ~isempty(ref_I) && isempty(outputSet)
                
                lb = zeros(obj.nI_fb);
                ub = zeros(obj.nI_fb);
                
                I2 = Polyhedron('lb', lb, 'ub', ub);
                I = Conversion.concatenatePolyhedron([ref_I I2]);
                
            end
            
            if ~isempty(ref_I) && ~isempty(outputSet)
               
                l = length(outputSet);
                nA = size(outputSet(1).A,2);
                I2 = [];
                for i=1:n

                    if l - obj.feedbackMap(i) <= 0
                        I2 = [I2 Polyhedron('lb', zeros(nA, 1), 'ub', zeros(nA, 1))];
                    else

                        I2 = [I2 outputSet(l - obj.feedbackMap(i))];

                    end                

                end

                I = Conversion.concatenatePolyhedron([ref_I I2]);
            end
            
            
            if isempty(ref_I) && ~isempty(outputSet)
                
                l = length(outputSet);
                nA = size(outputSet(1).A,2);
                I2 = [];
                for i=1:n

                    if l - obj.feedbackMap(i) <= 0
                        I2 = [I2 Polyhedron('lb', zeros(nA, 1), 'ub', zeros(nA, 1))];
                    else

                        I2 = [I2 outputSet(l - obj.feedbackMap(i))];

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

