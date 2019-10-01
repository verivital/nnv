classdef DLinearNNCS < handle
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
        % disturbance  --->|
        %              --->| plant ---> y(t) ---sampling--->y(k) 
        %             |                                       |
        %             |                                       |
        %             u(k) <---- controller |<---- v(k)-------|--- (reference inputs to the controller) 
        %                                   |<----- y(k)------|(output feedback)                                    
        
        
        % the input to neural net controller is grouped into 2 group
        % the first group contains all the reference inputs
           
        % the first layer weight matrix of the controller is decomposed into two
        % submatrices: W = [W1 W2] where
        %              W1 is conresponding to I1 = v[k] (the reference input)
        %              W2 is conresponding to I2 = y[k] (the feedback inputs)  
        
        % the reach set of the first layer of the controller is: 
        %              R = f(W1 * I1 + W2 * I2 + b), b is the bias vector of
        %              the first layer, f is the activation function
        
        nO = 0; % number of output
        nI = 0; % number of inputs = size(I1, 1) + size(I2, 1) to the controller
        nI_ref = 0; % number of reference inputs to the controller
        nI_fb = 0; % number of feedback inputs to the controller
        nI_plant_db = 0; % number of disturbance inputs to the plant 
             
        
        % for reachability analysis
        method = 'exact-star'; % by default
        plantReachSet = {};
        controllerReachSet = {};
        numCores = 1; % default setting, using single core for computation
        ref_I = []; % reference input set
        init_set = []; % initial set for the plant
        plant_db_set = []; % plant disturbance input set
        reachTime = 0; 
        
    end
    
    methods
        
        %constructor
        function obj = DLinearNNCS(controller, plant)
            % @controller: a neural net controller
            % @plant: a plant model (a LinearODE, DLinearODE or Neural net)
                       
            % author: Dung Tran
            % date: 9/30/2019
            
            
            if  ~isa(controller, 'FFNNS')
                error('The controller is not a feedforward neural network');
            end
            
            if ~controller.isPieceWiseNetwork
                error('The controller is not a piecewise network, i.e., there exists a layer whose activation function is not piecewise linear');
            end
            
            if ~isa(plant, 'DLinearODE')
                error('The plant is not a discrete linear system');
            end            
                        
            if plant.nO > controller.nI
                error('Inconsistency between number of feedback outputs and number of controller inputs');
            end
                        
            obj.controller = controller;
            obj.plant = plant;
            obj.nO = plant.nO; % number of outputs of the system == number of the plant's outputs
            obj.nI = controller.nI; % number of input to the controller
            obj.nI_fb = plant.nO; % number of feedback inputs to the controller
            obj.nI_ref = controller.nI - obj.nI_fb; % number of reference inputs to the controller
            obj.nI_plant_db = plant.nI - controller.nO; % number of disturbance inputs to the plant
            
        end
        
        
        
        % reach
        function [R, reachTime] = reach(obj, init_set, ref_input, method, numOfSteps, numCores)
            % @init_set: initial set of state, a star set
            % @ref_input: reference input, may be a vector or a star set
            % @method: 'exact-star' or 'approx-star'
            % @numOfSteps: number of steps
            % @numCores: number of cores used in computation
            
            % author: Dung Tran
            % date: 10/1/2019
            
            t = tic; 
            
            if ~isa(init_set, 'Star')
                error('Initial set is not a star set');
            end
            
            if numCores < 1
                error('Invalid number of cores used in computation');
            end
            
            if ~strcmp(method, 'exact-star') && ~strcmp(method, 'approx-star')
                error('Unknown reachability method, NNV currently supports exact-star and approx-star methods');
            end
            
            if numOfSteps < 1
                error('Invalid number of steps');
            end
            
            
            if ~isa(ref_input, 'Star') && size(ref_input, 2) ~= 1 && size(ref_input, 1) ~= obj.nI_ref
                error('Invalid reference input vector');
            end
            
            obj.ref_input = ref_input;
            obj.numCores = numCores;
            obj.method = method; 
            obj.init_set = init_set;
            
            obj.plantReachSet = cell(numOfSteps, 1);
            obj.controllerReachSet = cell(numOfSteps, 1);
            
            for k=1:numOfSteps
                
                X = obj.plantStepReach(k);
                
                if strcmp(obj.method, 'exact-star')
                
                    obj.plantReachSet{k} = X;
                
                elseif strcmp(obj.method, 'approx-star')
                    
                    B = Star.get_hypercube_hull(X);
                    obj.plantReachSet{k} = B.toStar; % combine all stars into a single over-approximate star
                    
                end
                
                obj.controllerReachSet{k} = obj.controllerStepReach(k);
                
                
            end
            
            reachTime = toc(t);
            
            obj.reachTime = reachTime; 
            
            R = obj.plantReachSet;
            
            
        end
        
        
        
        
        % plant step reach for step k
        function X = plantStepReach(obj, k)
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
                
                X = [];
                
                % compute exact reachable set of a plant using star set
                % X[k + 1] = A * X[k] + B * u[k]
                
                for i=1:nX
                    
                    V_X0 = obj.plant.A * X0(i).V;
                    
                    U0_i = U0{i}; % control set corresponding to the initial set of state X0(i)
                    
                    mU = length(U0_i); % number of control sets corresponding to the initial set of state X0(i)
                    
                    for j=1:mU
                        
                        V_U0 = obj.plant.B * U0_i(j).V;
                        
                        X_ij = Star(V_X0 + V_U0, U0_i(j).C, U0_i(j).d, U0_i(j).predicate_lb, U0_i(j).predicate_ub);
                        
                        X = [X X_ij];
                        
                    end
                    
                end
                                
            end
        end
        
        
        % get the input set to the plant
        function U = getPlantInputSet(obj, k)
            % @k: step
            % @U: control set to the plant
            
            % author: Dung Tran
            % date: 10/1/2019
            
            
            U1 = obj.controllerReachSet{k-1};
            
            if obj.nI_plant_db == 0
                U = U1;
            else
                
                n = length(U1);
                
                for i=1:n
                    
                    if ~isempty(obj.plant_db_set) && ~isa(obj.plant_db_set, 'Star')
                    
                        
                    
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


            I = obj.getControllerInputSet(k);
            n = length(I);
            for i=1:n
                U{i} = obj.controller.reach(I(i), 'exact-star', obj.numCores);
            end


        end
        
        % get input set for the controller at step k
        function I = getControllerInputSet(obj, k)
            % @k: step
            % @I: input set to the controller at step k
            
            % author: Dung Tran
            % date: 10/1/2019
            
            
            X = obj.plantReachSet{k};
            n = length(X);
            Y = [];
            for i=1:n
                Y = [Y X(i).affineMap(obj.plant.C)];
            end
            
            if obj.nI_ref == 0
                I = Y;
            else
                
                I = [];
                
                for i=1:n
                    
                    if isempty(obj.ref_I) % empty reference input is equal to zero vector
                    
                    I1 = Y(i).concatenate_with_vector(zeros(obj.nI_ref, 1));
                    
                    else 

                         if ~isa(obj.ref_I, 'Star') % referece input is just a vector
                             I1 = Y(i).concatenate_with_vector(obj.ref_I);
                         else
                             
                             if strcmp(obj.method, 'exact-star')
                                 
                                 if k == 1
                                    I1 = obj.ref_I.concatenate(Y(i)); % only do concatenate at the first step
                                 else                           
                                    % preseve relationship from second step 
                                    V1 = obj.ref_I.V;
                                    Z = zeros(obj.nI_ref, Y(i).nVar - obj.ref_I.nVar);
                                    V2 = Y(i).V;
                                    new_V = [V1 Z; V2];
                                    I1 = Star(new_V, Y(i).C, Y(i).d, Y(i).predicate_lb, Y(i).predicate_ub);

                                 end

                             elseif strcmp(obj.method, 'approx-star')
                                 
                                    I1 = obj.ref_I.concatenate(Y(i));
                                 
                             end
                             

                         end

                    end
 
                    I = [I I1];
                    
                end  
                
            end
            
            
        end
                
                
           
            
            
    end
        
    
    
    
    
end

