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
        
        % for reachability analysis
        method = 'exact-star'; % by default
        plantReachSet = {};
        controllerReachSet = {};
        numCores = 1; % default setting, using single core for computation
        ref_I = []; % reference input set
        init_set = []; % initial set for the plant
        reachTime = 0; 
        
    end
    
    methods % CONTRUCTOR AND REACH METHOD
        
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
            
            if plant.nI ~= controller.nO
                error('Inconsistency between the number of plant inputs and the number of controller outputs');
            end
                        
            obj.controller = controller;
            obj.plant = plant;
            obj.nO = plant.nO; % number of outputs of the system == number of the plant's outputs
            obj.nI = controller.nI; % number of input to the controller
            obj.nI_fb = plant.nO; % number of feedback inputs to the controller
            obj.nI_ref = controller.nI - obj.nI_fb; % number of reference inputs to the controller
            
        end
        
        
        
        % reach
        function [R, reachTime] = reach(varargin)
            % @init_set: initial set of state, a star set
            % @ref_input: reference input, may be a vector or a star set
            % @numOfSteps: number of steps
            % @method: 'exact-star' or 'approx-star'
            % @numCores: number of cores used in computation
            
            % author: Dung Tran
            % date: 10/1/2019
            
            
            switch nargin
                
                case 4
                    
                    obj = varargin{1};
                    init_set1 = varargin{2};
                    ref_input1 = varargin{3};
                    numOfSteps = varargin{4};
                    method1 = 'exact-star';
                    numCores1 = 1; 
                    
                case 5
                    obj = varargin{1};
                    init_set1 = varargin{2};
                    ref_input1 = varargin{3};
                    numOfSteps = varargin{4};
                    method1 = varargin{5};
                    numCores1 = 1;
                    
                case 6
                    obj = varargin{1};
                    init_set1 = varargin{2};
                    ref_input1 = varargin{3};
                    numOfSteps = varargin{4};
                    method1 = varargin{5};
                    numCores1 = varargin{6};
                    
                otherwise 
                    error('Invalid number of inputs, should be 3, 4, or 5');
            end
                
            
            
            t = tic; 
            
            if ~isa(init_set1, 'Star')
                error('Initial set is not a star set');
            end
            
            if numCores1 < 1
                error('Invalid number of cores used in computation');
            end
            
            if ~strcmp(method1, 'exact-star') && ~strcmp(method1, 'approx-star')
                error('Unknown reachability method, NNV currently supports exact-star and approx-star methods');
            end
            
            if numOfSteps < 1
                error('Invalid number of steps');
            end
            
            
            if ~isempty(ref_input1) && ~isa(ref_input1, 'Star') && size(ref_input1, 2) ~= 1 && size(ref_input1, 1) ~= obj.nI_ref
                error('Invalid reference input vector');
            end
            
            obj.ref_I = ref_input1;
            obj.numCores = numCores1;
            obj.method = method1; 
            obj.init_set = init_set1;
            
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
                Y = [Y X(i).affineMap(obj.plant.C, [])];
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
        
        
        
        % get output reach set
        
        function Y = getOutputReachSet(obj, map_mat, map_vec)
            % @map_mat: a mapping matrix
            % @map_vec: mapping vector
            % Y = map_mat * X + map_vec
            % @Y: a cell of output reach sets
            
            % author: Dung Tran
            % date: 10/2/2019
            
            
            if isempty(obj.plantReachSet )
                error('Plant reach set is empty, please perform reachability analysis first');
            end
            
            dim = size(obj.plant.A, 1); % dimension of the plant
            
            if size(map_mat, 2) ~= dim 
                error('Inconsistency between map_mat and dimension of the plant, map_mat should have %d columns', dim);
            end
            
            if size(map_mat, 1) > 3
                error('Plot only <= 3 dimensional outputs, the maximum allowable number of rows in map_mat is 3');
            end
            
            
            if ~isempty(map_vec) && (size(map_vec, 2) ~= 1)
                error('map vector should have one column');
            end
            
            if ~isempty(map_vec) && (size(map_vec, 1) ~= size(map_mat, 1))
                error('Inconsistent dimensions between map matrix and map vector');
            end
              
            % get output reach sets
          
            n = length(obj.plantReachSet);
            Y = cell(n, 1);
            for i=1:n
                
                Xi = obj.plantReachSet{i};
                m = length(Xi);
                Y1 = [];
                for j=1:m
                    Xij = Xi(j);
                    Y1 = [Y1 Xij.affineMap(map_mat, map_vec)];
                end
                
                Y{i} = Y1;
                
            end
            
        end
                     
            
    end
        
    
    methods % PLOT REACHABLE SETS
        
        
        function plotPlantReachSets(varargin)
            % @color: color
            % @x_id: id of x state
            % @y_id: id of y state
            % @z_id: id of z state
            
            % if plot in 2D, set one of three ids = []
            % if plot in 3D, three ids need to specified
            % if plot in 1D, set two of three ids = []
            % plot 1D: plot the ranges vs time steps
            
            % author: Dung Tran
            % date: 10/2/2019
            
                       
            switch nargin
                
                case 3
                    obj = varargin{1};
                    color = varargin{2};
                    dim = size(obj.plant.A,1); % dimension of the plant
                    x_id = varargin{3};
                    
                    if isempty(x_id) || x_id < 1 || x_id > dim
                        error('Invalid x index');
                    end
                    
                    option = 1; % plot range of the state with time steps 
                    
                case 4
                    obj = varargin{1};
                    color = varargin{2};
                    dim = size(obj.plant.A,1); % dimension of the plant
                    x_id = varargin{3};
                    y_id = varargin{4};
                    
                    if isempty(x_id) || x_id < 1 || x_id > dim
                        error('Invalid x index');
                    end
                    
                    if isempty(y_id) || y_id < 1 || y_id > dim
                        error('Invalid y index');
                    end
                    
                    if y_id == x_id 
                        error('x index and y index need to be different');
                    end
                    
                    option = 2; % plot 2D reachable set of X = (x_id, y_id)
                    
                case 5
                    obj = varargin{1};
                    color = varargin{2};
                    dim = size(obj.plant.A,1); % dimension of the plant
                    x_id = varargin{3};
                    y_id = varargin{4};
                    z_id = varargin{5};
                    
                    if isempty(x_id) || x_id < 1 || x_id > dim
                        error('Invalid x index');
                    end
                    
                    if isempty(y_id) || y_id < 1 || y_id > dim
                        error('Invalid y index');
                    end
                    
                    if isempty(z_id) || z_id < 1 || z_id > dim
                        error('Invalid z index');
                    end
                    
                    if y_id == x_id 
                        error('x index and y index need to be different');
                    end
                    
                    if y_id == z_id 
                        error('y index and z index need to be different');
                    end
                    
                    if z_id == x_id 
                        error('z index and x index need to be different');
                    end
                    
                    option = 3; % plot 3D reachable sets of X = (x_id, y_id, z_id)
            
                otherwise
                    error('Invalid number of inputs, should be 2, 3, or 4');
            end
            
            
            if isempty(obj.plantReachSet)
                error('Plant Reach Set is empty, please perform reachability analysis first.');
            end
            
            n = length(obj.plantReachSet); % number of steps
                        
            if option == 1 % plot ranges of specific state versus time steps
                
                X = obj.init_set;
                for i=1:n
                    X1 = Star.get_hypercube_hull(obj.plantReachSet{i});
                    X1 = X1.toStar;
                    X = [X X1];
                end
                T = 0:1:n;
                Star.plotRanges_2D(X, x_id, T, color);
                ax = gca; 
                ax.XTick = T;
                
            end
            
            if option == 2 % plot 2D reach set               
                
                X = obj.init_set; 
                for i=1:n
                    X = [X obj.plantReachSet{i}];
                end
                
                map = zeros(2, dim);
                map(1, x_id) = 1;
                map(2, y_id) = 1;
                                
                Y = [];
                for i=1:n+1
                    Y = [Y X(i).affineMap(map, [])]; 
                end
                
                Star.plots(Y, color);
                
            end
            
            if option == 3 % plot 3D reach set               
                
                X = obj.init_set; 
                for i=1:n
                    X = [X obj.plantReachSet{i}];
                end
                
                map = zeros(3, dim);
                map(1, x_id) = 1;
                map(2, y_id) = 1;
                map(3, z_id) = 1;
                
                Y = [];
                for i=1:n+1
                    Y = [Y X(i).affineMap(map, [])]; 
                end
                
                Star.plots(Y, color);
                
            end         
            
            
            
        end
        
        
        % plot controller reach sets
        function plotControllerReachSets(varargin)
            % @color: color
            % @x_id: id of first control input
            % @y_id: id of second control input
            % @z_id: id of third control input
            
            % if plot in 2D, set one of three ids = []
            % if plot in 3D, three ids need to specified
            % if plot in 1D, set two of three ids = []
            % plot 1D: plot the ranges vs time steps
            
            % author: Dung Tran
            % date: 10/2/2019
            
                       
            switch nargin
                
                case 3
                    obj = varargin{1};
                    color = varargin{2};
                    nU = size(obj.plant.B,2); % number of control inputs to the plant
                    x_id = varargin{3};
                    
                    if isempty(x_id) || x_id < 1 || x_id > nU
                        error('Invalid x index');
                    end
                    
                    option = 1; % plot range of the state with time steps 
                    
                case 4
                    obj = varargin{1};
                    color = varargin{2};
                    nU = size(obj.plant.B,2); % number of control inputs to the plant
                    x_id = varargin{3};
                    y_id = varargin{4};
                    
                    if isempty(x_id) || x_id < 1 || x_id > nU
                        error('Invalid x index');
                    end
                    
                    if isempty(y_id) || y_id < 1 || y_id > nU
                        error('Invalid y index');
                    end
                    
                    if y_id == x_id 
                        error('x index and y index need to be different');
                    end
                    
                    option = 2; % plot 2D reachable set of X = (x_id, y_id)
                    
                case 5
                    obj = varargin{1};
                    color = varargin{2};
                    nU = size(obj.plant.B,2); % number of control inputs to the plant
                    x_id = varargin{3};
                    y_id = varargin{4};
                    z_id = varargin{5};
                    
                    if isempty(x_id) || x_id < 1 || x_id > nU
                        error('Invalid x index');
                    end
                    
                    if isempty(y_id) || y_id < 1 || y_id > nU
                        error('Invalid y index');
                    end
                    
                    if isempty(z_id) || z_id < 1 || z_id > nU
                        error('Invalid z index');
                    end
                    
                    if y_id == x_id 
                        error('x index and y index need to be different');
                    end
                    
                    if y_id == z_id 
                        error('y index and z index need to be different');
                    end
                    
                    if z_id == x_id 
                        error('z index and x index need to be different');
                    end
                    
                    option = 3; % plot 3D reachable sets of X = (x_id, y_id, z_id)
            
                otherwise
                    error('Invalid number of inputs, should be 2, 3, or 4');
            end
            
            
            if isempty(obj.controllerReachSet)
                error('Controller Reach Set is empty, please perform reachability analysis first.');
            end
            
            n = length(obj.controllerReachSet); % number of steps
                        
            if option == 1 % plot ranges of specific state versus time steps
                
                U = [];
                for i=1:n
                    
                    U1 = Star.get_hypercube_hull(obj.controllerReachSet{i}{1});
                    U1 = U1.toStar;
                    U = [U U1];
                end
                T = 1:n;
                Star.plotRanges_2D(U, x_id, T, color);
                ax = gca; 
                ax.XTick = T;
                
            end
            
            if option == 2 % plot 2D reach set               
                
                U = []; 
                for i=1:n
                    U = [U obj.controllerReachSet{i}{1}];
                end
                
                map = zeros(2, dim);
                map(1, x_id) = 1;
                map(2, y_id) = 1;
                
                display(map);
                
                Y = [];
                n = length(U);
                for i=1:n
                    Y = [Y U(i).affineMap(map, [])]; 
                end
                
                Star.plots(Y, color);
                
            end
            
            if option == 3 % plot 3D reach set               
                
                U = [];
                for i=1:n
                    U = [U obj.controllerReachSet{i}{1}];
                end
                
                map = zeros(3, dim);
                map(1, x_id) = 1;
                map(2, y_id) = 1;
                map(3, z_id) = 1;
                
                Y = [];
                n = length(U);
                for i=1:n
                    Y = [Y U(i).affineMap(map, [])]; 
                end
                
                Star.plots(Y, color);
                
            end         
              
            
        end
        
        
        % plot output reach set
        % output reach set is derived by mapping the state reach set by
        % a maping matrix, M 
        function plotOutputReachSets(obj, color, map_mat, map_vec)
            % @map_mat: a mapping matrix
            % @map_vec: mapping vector
            % Y = map_mat * X + map_vec
            % @color: color
            
            % author: Dung Tran
            % date: 10/2/2019
            
            
            Y = obj.getOutputReachSet(map_mat, map_vec);
            
            n = length(Y);
            
            % plot output reach sets
            
            option = size(map_mat, 1);
            
            G = obj.init_set.affineMap(map_mat, map_vec);
            
            if option == 1 % plot 1D, output versus time steps        
                                
                for i=1:n
                    G1 = Star.get_hypercube_hull(Y{i});
                    G1 = G1.toStar;
                    G = [G G1];
                end
                T = 0:1:n;
                Star.plotRanges_2D(G, 1, T, color);
                ax = gca; 
                ax.XTick = T;
                
            end
            
            if option == 2 || option == 3 % plot 2D or 3D
                
                for i=1:n
                    G = [G Y{i}];
                end
                
                Star.plots(G, color);
                
            end  
            
            
        end
          
        
    end
    
    
    
    methods % VERIFICATION METHOD
        
        
        function [safe, counterExamples] = verify(obj, unsafe_mat, unsafe_vec)
            % @unsafe_mat: unsafe matrix
            % @unsafe_vec: unsafe vector            
            % Usafe region is defined by: x: unsafe_mat * x <= unsafe_vec
            
            % @safe: = 0: unsafe
            %        = 1: safe
            %        = 2: unknown (due to conservativeness)
            % @counterExamples: an array of star set counterExamples
            
            
            % author: Dung Tran
            % date: 10/2/2019
            
            
            if isempty(obj.plantReachSet)
                error('Plant reachable set is empty, please do reachability analysis first');
            end
            
            dim = size(obj.plant.A, 1); 
            
            if (size(unsafe_mat, 2) ~= dim) 
                error('Inconsistent dimensions between the unsafe matrix and the plant');
            end
            
            if size(unsafe_vec, 2) ~= 1
                error('Unsafe vector should have one column');
            end
            
            if size(unsafe_vec, 1) ~= size(unsafe_mat, 1)
                error('Inconsistent dimension between unsafe matrix and unsafe vector');
            end
            
            X = obj.plantReachSet; 
            n = length(X);
            % flatten reach sets
            Y = [];
            for i=1:n
                Y = [Y X{i}];
            end
            
            % check safety
            m = length(Y);
            G = [];
            for i=1:m
                G1 = Y(i).intersectHalfSpace(unsafe_mat, unsafe_vec);
                if ~G1.isEmptySet
                    G = [G G1];
                end
            end
            
            if isempty(G)
                safe = 1;
            else
                
                if strcmp(obj.method, 'exact-star')
                    safe = 0;
                    % construct a set of counter examples                    
                    n = length(G);
                    counterExamples = [];
                    for i=1:n
                        C1 = Star(obj.init_set.V, G(i).C, G(i).d, G(i).predicate_lb, G(i).predicate_ub);
                        counterExamples = [counterExamples C1];
                    end
                    
                    
                elseif strcmp(obj.method, 'approx-star')
                    safe = 2; % unknown due to over-approximation error, we can try using simulation to falsify the property
                    counterExamples = [];
                end
            
            
            end
            
            
            if safe == 0
                fprintf('\n\nThe neural network control system is unsafe, NNV produces %d star sets of counter-examples', n);
            elseif safe == 1
                fprintf('\n\nThe neural network control system is safe');
            elseif safe == 2
                fprintf('\n\n The safety of the neural network control system is unknown');
                fprintf('\nYou can try falsification method using simulation to find counter-examples');
            end
            
            
        end
        
        
        
    end
    
    
    
    
end

