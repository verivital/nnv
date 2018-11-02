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
        
        % feedbackMap = [1;2], a 2 x 1 matrix, means that:
        % the output feedback to the controller are: [y[k-1]; y[k-2]]        
        
        % the first layer weight matrix of the controller is decomposed into two
        % submatrices: W = [W1 W2] where
        %              W1 is conresponding to I1 = v[k] (the reference input)
        %              W2 is conresponding to I2 = [y[k-1]; y[k-2]] (the feedback inputs)  
        
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
     function P = reachPolyhedron(obj, init_set, ref_inputSet, N)
         % @init_set: the initial set of condition for the plant
         % @ref_inputSet: the reference input set applied to the controller
         % @N: number of steps 
         % @R: the reachable set of the plant
         %     we get the output reachable set by mapping R on the
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
             P = Rx;
         end         
         
         if N > 1
             
             % first step
             Rx = obj.plant.stepReachPolyhedron(init_set, []); % state reachable set
             P = Rx;
             
             
             for i=2:N
                 % next input set for the controller
                 
                 
                 
             end
             
             
         end
         
         
         
         
     end
     
     
        
        
    end
    
    
    methods(Static)
        
         % step Reach controller
         function U = stepReachController(controller, ref_input, fb_input)
             % @controller: the controller neural network
             % @ref_input: reference input set
             % @fb_input: feedback input set
            
             

         end

        
                           
    end
    
end

