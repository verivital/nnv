classdef LinearNNCS < handle
    %Class for Linear Neural Network Control Systems
    %   Dung Tran
    %   Date: 9/30/2019
    
    % TODO: need support feedback y(k), y(k-d)
    
    properties
        controller = []; % a feedfoward neural network
        plant = []; % linear plant model, DLinearODE or LinearODE
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
        
        plantReachSet = {};
        controllerReachSet = {};
        
        numCores = 1; % default setting, using single core for computation
        
    end
    
    methods
        
        %constructor
        function obj = LinearNNCS(controller, plant, feedbackMap)
            % @controller: a neural net controller
            % @plant: a plant model (a LinearODE, DLinearODE or Neural net)
            % @feedbackMap: a feedback map from outputs of the plant to the
            % input of the controller
            
            % author: Dung Tran
            % date: 9/30/2019
            
            
            if  ~isa(controller, 'FFNNS')
                error('The controller is not a feedforward neural network');
            end
            
            if ~controller.isPieceWiseNetwork
                error('The controller is not a piecewise network, i.e., there exists a layer whose activation function is not piecewise linear');
            end
            
            if ~isa(plant, 'LinearODE') && ~isa(plant, 'DLinearODE')
                error('The plant is not a linear ode system in discrete or continuous time');
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
        
        
        
        
        % plant step reach for step k
        function X = plantStepReachDiscrete(obj, k)
            % @k: step
            % @X: plant reach set
            
            % author: Dung Tran
            % date: 9/30/2019
            
            
            
            if k==1
                X = obj.init_set.affineMap(obj.plant.A, []);
            end
            
            if k>1
                X0 = obj.plantReachSet{k-1};
                U0 = obj.controllerReachSet{k-1};
                
                nX = length(X0);
                [nU, mU] = size(U0);
                
                if nX ~= nU
                    error('Inconsistent between initial set of state and control set');
                end
                
                X = [];
                
                % compute exact reachable set of a plant using star set
                % X[k + 1] = A * X[k] + B * u[k]
                
                for i=1:nX
                    
                    V_X0 = obj.plant.A * X0(i).V; 
                    
                    for j=1:mU
                        
                        V_U0 = obj.plant.B * U0{i, j}.V;
                        
                        X_ij = Star(V_X0 + V_U0, U0{i, j}.C, U0{i, j}.d, U0{i, j}.predicate_lb, U0{i, j}.predicate_ub);
                        
                        X = [X X_ij];
                        
                    end
                    
                end
                                
            end
        end
            
            
        % controller step Reach for step k
        function U = controllerStepReach(obj, k)
            % @k: step
            % @U: controller reach set

            % author: Dung Tran
            % date: 10/1/2019


            X = obj.plantReachSet{k};
            n = length(X);
            U = cell(n,1);

            for i=1:n
                U{i} = obj.controller.reach(X(i), 'exact-star', obj.numCores);
            end


        end
        
        % get input set for the controller at step k
        function I = getControllerInputSet(obj, k)
            % @k: step
            % @I: input set to the controller at step k
            
            % author: Dung Tran
            % date: 10/1/2019
            
            
            
            
            
        end
                
                
           
            
            
    end
        
    
    
    
    
end

