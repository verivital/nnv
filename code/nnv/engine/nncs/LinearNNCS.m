classdef LinearNNCS < handle
    %Class for Linear Neural Network Control Systems
    %   Dung Tran
    %   Date: 9/30/2019
    
    % TODO: need support feedback y(k), y(k-d)
    
    properties
        controller = []; % a feedfoward neural network
        plant = []; % linear plant model
        % nerual network control system architecture
        %
        % nerual network control system architecture
        %
        %              ---> plant ---> y(t) ---sampling--->y(k) 
        %             |                                       |
        %             |                                       |
        %             u(k) <---- controller |<---- y(k)-----(output feedback) 
        %                                   |<----- v(k)------(reference input)   
        %                                   
        
        nO = 0; % number of output
        nI = 0; % number of inputs = size(I1, 1) + size(I2, 1) to the controller
        nI_ref = 0; % number of reference inputs to the controller
        nI_fb = 0; % number of feedback inputs to the controller
        
        % for reachability analysis
        method = 'exact-star'; % by default
        plantReachMethod = 'direct'; % by default, direct method, compute e^Ah
        %                = 'ode45'  : using Ode45 solver
        % Look at LinearODE class for the detail
        
        transPlant = [];  % a transformed plant for reachability analysis  
        % A_new = [A B; 0 0]; C_new = [C 0]
        % transPlant: dot{x_new} = A_new * x_new; y_new = C_new * x_new
        
        plantReachSet = {};
        plantIntermediateReachSet = {};
        plantNumOfSimSteps = 20; % default number of simulation steps for the plant between two control step k-1 and k
        controlPeriod = 0.1; % default control period
        controllerReachSet = {};
        numCores = 1; % default setting, using single core for computation
        ref_I = []; % reference input set
        init_set = []; % initial set for the plant
        reachTime = 0;
        
        % for simulation
        simTraces = {}; % simulation trace
        controlTraces = {}; % control trace
        
        % use for falsification
        falsifyTraces = {};
        falsifyTime = 0;
        
    end
    
    methods % CONTRUCTOR AND REACH METHOD
        
        %constructor
        function obj = LinearNNCS(controller, plant)
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
            
            if ~isa(plant, 'LinearODE')
                error('The plant is not a linear system');
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
            obj.plantNumOfSimSteps = obj.plant.numReachSteps;
            obj.controlPeriod = obj.plant.controlPeriod;
            new_A = [obj.plant.A obj.plant.B; zeros(obj.plant.nI, obj.plant.dim) zeros(obj.plant.nI, obj.plant.nI)];
            new_C = [eye(obj.plant.dim) zeros(obj.plant.dim, obj.plant.nI)];
            obj.transPlant = LinearODE(new_A, [], new_C, [], obj.plant.controlPeriod, obj.plant.numReachSteps);
            
        end
        
        
        % start parallel pool for computing
        function start_pool(obj)
            
            if obj.numCores > 1
                poolobj = gcp('nocreate'); % If no pool, do not create new one.
                if isempty(poolobj)
                    parpool('local', obj.numCores); 
                else
                    if poolobj.NumWorkers ~= obj.numCores
                        delete(poolobj); % delete the old poolobj
                        parpool('local', obj.numCores); % start the new one with new number of cores
                    end                    
                end
            end   
            
        end
        
        
        % reach
        function [R, reachTime] = reach(obj, reachPRM)
            % @reachPRM: reachability parameters containing following
            % inputs:
            %       1) @init_set: initial set of state, a star set
            %       2) @ref_input: reference input, may be a vector or a star set
            %       3) @numOfSteps: number of steps
            %       4) @method: 'exact-star' or 'approx-star'
            %       5) @numCores: number of cores used in computation
            %       6) @numOfSimSteps: number of sim steps to compute reachable set
            %          for the plant
            % NOTE***: parallel computing may not help due to
            % comuninication overhead in the computation
            
            
            % author: Dung Tran
            % date: 10/1/2019
            % update: 3/15/2020
            
            init_set1 = reachPRM.init_set;
            ref_input1 = reachPRM.ref_input;
            numOfSteps = reachPRM.numSteps; 
            method1 = reachPRM.reachMethod;
            numCores1 = reachPRM.numCores;
                       
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
            
            obj.plantReachSet = cell(1, numOfSteps);
            obj.plantIntermediateReachSet = cell(1,numOfSteps);
            obj.controllerReachSet = cell(1, numOfSteps);
                            
            if obj.numCores > 1
                obj.start_pool;
            end
            
            t = tic; 
            
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
        
        % live reachability analysis, plot reachable set on the fly
        % produce video for the analysis
        function [R, reachTime] = reachLive(varargin)
            % @init_set: initial set of state, a star set
            % @ref_input: reference input, may be a vector or a star set
            % @numOfSteps: number of steps
            % @method: 'exact-star' or 'approx-star'
            % @numCores: number of cores used in computation
            % NOTE***: parallel computing may not help due to
            % comuninication overhead in the computation
            % @map_mat: output mapping matrix 
            % @map_vec: output mapping bias vector
            % *** We plot y = map_mat * x + map_vec on-the-fly
            % @color: color for plotting
            
            % author: Dung Tran
            % date: 10/1/2019
            % update: 11/6/2019
            
            obj = varargin{1};
            option.initSet = [];
            option.refInput = [];
            option.reachMethod = 'approx-star';
            option.numCores = 1;
            option.plantNumOfSimSteps = 20;
            map_mat = zeros(1, obj.plant.dim);
            map_mat(1) = 1;
            option.outputMatrix = map_mat;
            option.outputVector = [];
            option.outputSetColor = 'blue';
            option.boundaryMatrix = [];
            option.boundaryVector = [];
            option.boundarySetColor = 'red'; 
            option.figureTitle = 'Reachable Sets';
            option.figureXLabel = '';
            option.figureYLabel = '';
            option.videoRecord = true;
            option.videoName = 'reachVideo';
            option.videoFrameRate = 4; % frame per second
            
            n = length(varargin);
            if n < 4
                error('Not enough inputs for analysis');
            end
            
            option.initSet = varargin{2};
            if ~isa(option.initSet, 'Star')
                error('Initial set is not a star set');
            end
            
            option.refInput = varargin{3};
            if ~isempty(option.refInput) && ~isa(option.refInput, 'Star') && size(option.refInput, 2) ~= 1 && size(option.refInput, 1) ~= obj.nI_ref
                error('Invalid reference input vector');
            end
            
            option.numOfSteps = varargin{4};
            if option.numOfSteps < 1
                error('Invalid number of control steps');
            end
 
            for i=5:n
                if ischar(varargin{i})
                    
                    % parse reach method
                    if strcmp(varargin{i}, 'reachMethod')                       
                        if ~strcmp(varargin{i+1}, 'approx-star') || ~strcmp(varargin{i+1}, 'exact-star')
                            error('Unknown reachability method');
                        else
                            option.reachMethod = varargin{i+1};
                        end
                    end
                                        
                    % parse numCores
                    if strcmp(varargin{i}, 'numCores')
                        if varargin{i+1} < 1 
                            error('Invalid number of cores used for computation');
                        else
                            option.numCores = varargin{i+1};
                        end
                    end
                    
                    % parse plantNumOfSimSteps
                    if strcmp(varargin{i}, 'plantNumOfSimSteps')
                        if varargin{i+1} < 1 
                            error('Invalid number of simulation steps for reachability of the plant');
                        else
                            option.plantNumOfSimSteps = varargin{i+1};
                        end
                    end
                    
                    % parse output matrix
                    if strcmp(varargin{i}, 'outputMatrix')
                        if ~ismatrix(varargin{i+1}) 
                            error('Invalid output Matrix');
                        else 
                            option.outputMatrix = varargin{i+1};
                        end
                    end
                    
                    % parse output vector
                    if strcmp(varargin{i}, 'outputVector')
                        if ~ismatrix(varargin{i+1}) 
                            error('Invalid output Vector');
                        else 
                            option.outputVector = varargin{i+1};
                        end
                    end
                    
                    % parse boundary matrix
                    if strcmp(varargin{i}, 'boundaryMatrix')
                        if ~ismatrix(varargin{i+1}) 
                            error('Invalid boundary Matrix');
                        else 
                            option.boundaryMatrix = varargin{i+1};
                        end
                    end
                    
                    % parse boundary matrix
                    if strcmp(varargin{i}, 'boundaryVector')
                        if ~ismatrix(varargin{i+1}) 
                            error('Invalid boundary Vector');
                        else 
                            option.boundaryVector = varargin{i+1};
                        end
                    end
                    
                    % parse outputReachSet color
                    if strcmp(varargin{i}, 'outputSetColor')
                        option.outputSetColor = varargin{i+1};
                    end
                    
                    % parse boundaryReachSet color
                    if strcmp(varargin{i}, 'boudarySetColor')
                        option.boundarySetColor = varargin{i+1};
                    end
                    
                    % parse figureTitle
                    if strcmp(varargin{i}, 'figureTitle')
                        option.figureTitle = varargin{i+1};
                    end
                    
                     % parse figureXLabel
                    if strcmp(varargin{i}, 'figureXLabel')
                        option.figureXLabel = varargin{i+1};
                    end
                    
                    
                    % parse figureYLabel
                    if strcmp(varargin{i}, 'figureYLabel')
                        option.figureYLabel = varargin{i+1};
                    end
                    
                    % parse videoRecord
                    if strcmp(varargin{i}, 'videoRecord')
                        option.videoRecord = varargin{i+1};
                    end
                    
                    % parse videoName
                    if strcmp(varargin{i}, 'videoName')
                        option.videoName = varargin{i+1};
                    end
                    
                    % parse videoFrameRate
                    if strcmp(varargin{i}, 'videoFrameRate')
                        option.videoFrameRate = varargin{i+1};
                    end
                    
                end

            end
                      
            obj.ref_I = option.refInput;
            obj.numCores = option.numCores;
            obj.method = option.reachMethod; 
            obj.init_set = option.initSet;
            
            obj.plantNumOfSimSteps = option.plantNumOfSimSteps;
            obj.plantReachSet = cell(option.plantNumOfSimSteps, 1);
            obj.controllerReachSet = cell(option.plantNumOfSimSteps, 1);
            
                
            if obj.numCores > 1
                obj.start_pool;
            end
            
            
            if option.videoRecord                
                option.reachVideo = VideoWriter(option.videoName); % Name it.
                option.reachVideo.FrameRate = option.videoFrameRate; % How many frames per second.
                open(option.reachVideo);
            end

            t = tic; 
            
            for k=1:option.numOfSteps
                
                X = obj.plantStepReach(k);
                
                if strcmp(obj.method, 'exact-star')
                
                    obj.plantReachSet{k} = X;
                
                elseif strcmp(obj.method, 'approx-star')
                    
                    B = Star.get_hypercube_hull(X);
                    obj.plantReachSet{k} = B.toStar; % combine all stars into a single over-approximate star
                    
                end
                
                obj.controllerReachSet{k} = obj.controllerStepReach(k); 
                
                % live plotting reachable set
                
                if size(map_mat, 1) == 1
                    
                    obj.plotStepOutputReachSets(option, k);
                    
                    N = fix(k/10);
                   
                    if N == 0
                        times = 0: obj.controlPeriod: k*obj.controlPeriod;
                    else
                        times = 0: N*obj.controlPeriod: k*obj.controlPeriod;
                    end
                    ax = gca;
                    ax.XTick = times;
                    hold on;      
                                        
                elseif (size(map_mat, 1) == 2) || (size(map_mat, 1) == 3)
                    
                    obj.plotStepOutputReachSets(option, k);
                    
                elseif size(indexes, 1) > 3
                    error('NNV plots only three-dimensional output reachable set, please limit number of rows in map_mat <= 3');
                end
                
                 % make video
                if option.videoRecord
                    frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
                    writeVideo(option.reachVideo, frame);
                end   
                

            end
            
            if option.videoRecord
                close(option.reachVideo);
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
            
            
            
            h = obj.controlPeriod/obj.plantNumOfSimSteps; % simulation timestep for the plant model 
            if k==1
                
                % step 1: construct initial set for the transformed plant
                new_V = [obj.init_set.V; zeros(obj.plant.nI, obj.init_set.nVar + 1)];
                trans_init_set = Star(new_V, obj.init_set.C, obj.init_set.d, obj.init_set.predicate_lb, obj.init_set.predicate_ub);
                % step 2: perform the reachability analysis for the time elapse in the first step    
                
                 X_imd_trans = obj.transPlant.simReach(obj.plantReachMethod, trans_init_set, [], h, obj.plantNumOfSimSteps, []);
                 X_imd = [];
                 for i=1:obj.plantNumOfSimSteps +1
                     X_imd = [X_imd X_imd_trans(i).affineMap(obj.transPlant.C, [])];
                 end
                 obj.plantIntermediateReachSet{k} = cell(1,1);
                 obj.plantIntermediateReachSet{k}{1,1} = {X_imd};
                 X = X_imd(obj.plantNumOfSimSteps+1);
            end
            
            if k>1
                
                X0 = obj.plantReachSet{k-1};
                U0 = obj.controllerReachSet{k-1};
                
                nX = length(X0);
                                                                
                % compute exact reachable set of the continuous plant from t[k] to t[k+1] using star set
                
                X = [];
                obj.plantIntermediateReachSet{k} = cell(1, nX); 
                for i=1:nX                                        
                    U0_i = U0{i}; % control set corresponding to the initial set of state X0(i)
                    mU = length(U0_i); % number of control sets corresponding to the initial set of state X0(i)
                    obj.plantIntermediateReachSet{k}{i} = cell(mU, 1);
                    for j=1:mU
                        new_V = [X0(i).V; U0_i(j).V];
                        trans_init_set = Star(new_V, U0_i(j).C, U0_i(j).d, U0_i(j).predicate_lb, U0_i(j).predicate_ub);
                        X_imd = obj.transPlant.simReach(obj.plantReachMethod, trans_init_set, [], h, obj.plantNumOfSimSteps, []);
                        X1 = [];
                        for l=1:obj.plantNumOfSimSteps + 1
                            X1 = [X1 X_imd(l).affineMap(obj.transPlant.C, [])]; 
                        end
                        X = [X X1(obj.plantNumOfSimSteps+1)]; % store plant reach set at t[k], which is the initial set for the next control period t[k+1]
                        obj.plantIntermediateReachSet{k}{i}{j} = X1; % store intermediate plant reach set
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
            
            
            if isempty(obj.plantIntermediateReachSet )
                error('Plant reach set is empty, please perform reachability analysis first');
            end
            
            dim = obj.plant.dim; % dimension of the plant
            
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
            n = length(obj.plantIntermediateReachSet);
            
            Y = cell(1, n);
            for i=1:n
                Y{i} = obj.getStepOutputReachSet(map_mat, map_vec, i);
            end
            
        end
        
        % get step output reachable set
        function Y = getStepOutputReachSet(obj, map_mat, map_vec, k)
            % @map_mat: a mapping matrix
            % @map_vec: mapping vector
            % Y = map_mat * X + map_vec
            % @Y: a cell of output reach sets
            % @k: step index
            
            % author: Dung Tran
            % date: 10/2/2019
            
            
            if isempty(obj.plantIntermediateReachSet )
                error('Plant reach set is empty, please perform reachability analysis first');
            end
            
            dim = obj.plant.dim; % dimension of the plant
            
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
          
               
            X = obj.plantIntermediateReachSet{k};
            nX = length(X);
            
            Y = cell(1, obj.plantNumOfSimSteps+1);
            
            for l=1:obj.plantNumOfSimSteps+1
                Y1 = [];
                for i=1:nX
                    Xi = X{i};
                    nXi = length(Xi);
                    for j=1:nXi
                        Xij = Xi{j};
                        Yijl = Xij(l).affineMap(map_mat, map_vec);
                        Y1 = [Y1 Yijl];
                    end
                end
                Y{l} = Y1;
            end
            
        end
                     
            
    end
        
    
    methods % PLOT REACHABLE SETS
        
        
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
            
            n = length(Y); % number of control periods
            
            % plot output reach sets
            
            option = size(map_mat, 1);
            
            I0 = obj.init_set.affineMap(map_mat, map_vec);
            G = I0.getBox;
            h = obj.controlPeriod / obj.plantNumOfSimSteps;
                        
            if option == 1 % plot 1D, output versus time steps        
                
                for i=1:n
                    Y1 = Y{i}; 
                    G1 = [];
                    for j=1:obj.plantNumOfSimSteps+1
                        B = Star.get_hypercube_hull(Y1{j});
                        G1 = [G1 B.toStar];
                    end
                    
                    T1 = (i-1)*obj.controlPeriod:h:i*obj.controlPeriod;
                    Star.plotRanges_2D(G1, 1, T1, color);
                    hold on;
                end
                
                T = 0:obj.controlPeriod:n*obj.controlPeriod;
                ax = gca; 
                ax.XTick = T;
                
            end
            
            if option == 2 || option == 3 % plot 2D or 3D
                
                for i=1:n
                    Y1 = Y{i}; 
                    for j=1:obj.plantNumOfSimSteps+1
                        B = Star.get_hypercube_hull(Y1{j});
                        G = [G B];
                    end
                end
                
                if option == 2
                    Box.plotBoxes_2D_noFill(G, 1, 2, color);
                elseif option == 3
                    Box.plotBoxes_3D_noFill(G, 1, 2, 3, color);
                end
                
            end
            
            if option > 3
                error('We can plot only 3-dimensional output set, please limit the number of row of the mapping matrix to <= 3');
            end
            
            
        end
        
        
        % plot step output reachable sets
        function plotStepOutputReachSets(obj, option, k)
            % @out_mat: a mapping matrix
            % @out_vec: mapping vector
            % Y = out_mat * X + out_vec
            % @unsafe_mat: unsafe region matrix
            % @unsafe_vec: unsafe region vector
            % @color: color
            % @k: step index
            
            % author: Dung Tran
            % date: 10/2/2019
            % update: 11/2/2019
            
            
            Y = obj.getStepOutputReachSet(option.outputMatrix, option.outputVector, k); % output set at step k
            
            if ~isempty(option.boundaryMatrix) && ~isempty(option.boundaryVector)               
                US = obj.getStepOutputReachSet(option.boundaryMatrix, option.boundaryVector, k); % unsafe region output
            else
                US = [];
            end
            
            h = obj.controlPeriod / obj.plantNumOfSimSteps;
                       
            % plot output reach sets
            
            Dim = size(option.outputMatrix, 1);
            
                       
            if Dim == 1 % plot 1D, output versus time steps

               for i=2:obj.plantNumOfSimSteps + 1
                    B = Star.get_hypercube_hull(Y{i});
                    ymin = B.lb;
                    ymax = B.ub;
                    y = [ymin ymin ymax ymax ymin];
                    xmin = (k-1)*obj.controlPeriod + (i-2)*h;
                    xmax = xmin + h;
                    x = [xmin xmax xmax xmin xmin];
                    plot(x, y, option.outputSetColor);
                    hold on;
                    
                    if ~isempty(US)
                        B = Star.get_hypercube_hull(US{i});
                        ymin = B.lb;
                        ymax = B.ub;
                        y = [ymin ymin ymax ymax ymin];
                        xmin = (k-1)*obj.controlPeriod + (i-2)*h;
                        xmax = xmin + h;
                        x = [xmin xmax xmax xmin xmin];
                        plot(x, y, option.boundarySetColor);
                        
                    end
                    
                    ax = gca;
                    title(option.figureTitle);
                    xlabel(ax, option.figureXLabel);
                    ylabel(ax, option.figureYLabel);         
                    
                    pause(0.25);
                     % make video
                    if option.videoRecord
                        frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
                        writeVideo(option.reachVideo, frame);
                    end   
                    
                    
                    
               end                    

            end
            
            if Dim == 2 || Dim == 3 % plot 2D or 3D
                
                 for i=2:obj.plantNumOfSimSteps + 1
                     
                     B = Star.get_hypercube_hull(Y{i});
                     if Dim == 2
                        Star.plotBoxes_2D_noFill(B.toStar, 1, 2, option.outputSetColor);
                     else
                        Star.plotBoxes_3D_noFill(B.toStar, 1, 2, 3, option.outputSetColor);
                     end
                     hold on;
                    
                     if ~isempty(US)
                         B = Star.get_hypercube_hull(US{i});
                          if Dim == 2                             
                              Star.plotBoxes_2D_noFill(B.toStar, 1, 2, option.boundarySetColor);
                          else
                              Star.plotBoxes_3D_noFill(B.toStar, 1, 2, 3, option.boundarySetColor);
                          end
                         hold on;
                     end
                     
                     ax = gca;
                     title(option.figureTitle);
                     xlabel(ax, option.figureXLabel);
                     ylabel(ax, option.figureYLabel); 
                     
                     pause(0.25);
                     
                      % make video
                     if option.videoRecord
                         frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
                         writeVideo(option.reachVideo, frame);
                     end                
                     
                 end       
                
            end
            
            if Dim > 3
                error('We can plot only 3-dimensional output set, please limit the number of row of the mapping matrix to <= 3');
            end
            
            
            
        end
          
        
    end
    
    
    
    methods % VERIFICATION METHOD
        
        
        function [safe, counterExamples, verifyTime] = verify(obj, reachPRM, unsafeRegion)
            % @reachPRM: reachability parameters consists of following
            % inputs: 
            %       1) @reachPRM.init_set: initial set
            %       2) @reachPRM.ref_input: reference input
            %       3) @reachPRM.numSteps: number of steps
            %       4) @reachPRM.reachMethod: method for reachability analysis
            %       5) @reachPRM.numCores: number of cores used in reachability analysis
            %       6) @reachPRM.plantNumSimSteps: number of simulation
            %          steps for the plant in 1 control period
            %       6) @reachPRM.controlPeriod: controler period 
            % @unsafeRegion: a Halfpsace object
            % Usafe region is defined by: y: unsafe_mat * x <= unsafe_vec
            
            % @safe: = unsafe
            %        = safe
            %        = unknown (due to conservativeness)
            % @counterExamples: an array of star set counterExamples or
            %                   falsified input points
            % @verifyTime: verification time
            
            % author: Dung Tran
            % date: 10/2/2019
            % update: 3/15/2020
            
                        
            unsafe_mat = unsafeRegion.G;
            unsafe_vec = unsafeRegion.g;
            falsifyPRM.init_set = reachPRM.init_set;
            falsifyPRM.ref_input = reachPRM.ref_input;
            falsifyPRM.numSteps = reachPRM.numSteps;
            falsifyPRM.numTraces = 1000;
            falsifyPRM.controlPeriod = obj.plant.controlPeriod;
            falsifyPRM.unsafeRegion = unsafeRegion;
                        
            t = tic; 
                        
            dim = obj.plant.dim; 
            
            if (size(unsafe_mat, 2) ~= dim) 
                error('Inconsistent dimensions between the unsafe matrix and the plant');
            end
            
            if size(unsafe_vec, 2) ~= 1
                error('Unsafe vector should have one column');
            end
            
            if size(unsafe_vec, 1) ~= size(unsafe_mat, 1)
                error('Inconsistent dimension between unsafe matrix and unsafe vector');
            end
            
            obj.reach(reachPRM); % perform reachability analysis
            X = obj.plantReachSet; 
            n = length(X);
            % flatten reach sets
            Y = obj.init_set;
            for i=1:n
                Y = [Y X{i}];
            end
            
            % check safety
            m = length(Y);
            G = [];
            for i=1:m
                G1 = Y(i).intersectHalfSpace(unsafe_mat, unsafe_vec);
                if ~isempty(G1)
                    G = [G G1];
                end
            end
            
            if isempty(G)
                safe = 'SAFE';
                counterExamples = [];
            else
                
                if strcmp(obj.method, 'exact-star')
                    safe = 'UNSAFE';
                    % construct a set of counter examples                    
                    n = length(G);
                    counterExamples = [];
                    for i=1:n
                        C1 = Star(obj.init_set.V, G(i).C, G(i).d, G(i).predicate_lb, G(i).predicate_ub);
                        counterExamples = [counterExamples C1];
                    end
                    
                    
                elseif strcmp(obj.method, 'approx-star')
                    [safe, counterExamples, ~] = obj.falsify(falsifyPRM);
                end
            
            
            end
            
            
            if strcmp(safe, 'UNSAFE')
                fprintf('\n\nThe neural network control system is unsafe, NNV produces counter-examples');
            elseif strcmp(safe, 'SAFE')
                fprintf('\n\nThe neural network control system is safe');
            elseif strcmp(safe, 'UNKNOWN')
                fprintf('\n\n The safety of the neural network control system is unknown');
                fprintf('\nYou can try falsification method using simulation to find counter-examples');
            end
            
            verifyTime = toc(t);
            
            
        end
        
        
        
    end
    
    
    methods %SIMULATION
        
        
        % simulate the system
        function [simTrace, controlTrace, simTime] = simulate(obj, init_state, ref_input, numSteps)
           % @init_state: a vector of initial states
           % @ref_input: a vector of reference inputs
           % @numSteps: number of simulation steps
           % @simTrace: simulation trace
           % @controlTrace: control trace
           % @simTime: simulation time
           
           % author: Dung Tran
           % date: 10/2/2019
           
           t = tic;
           
           dim = obj.plant.dim; % dimension of the plant
      
           if size(init_state, 1) ~= dim
               error('Inconsistent between init_state vector and plant dimension');
           end
           
           if size(init_state, 2) ~= 1
               error('Invalid init_state vector, should have one column');
           end
           
           if ~isempty(ref_input) && (size(ref_input, 2) ~= 1)
               error('Invalid ref_input vector, should have one column');
           end
           
           if numSteps < 1
               error('Invalid number of steps, should be >= 1');
           end
           
           if ~isempty(ref_input) && (size(ref_input, 1) ~= obj.nI_ref)
               error('Invalid ref_input vector, should have %d elements', obj.nI_ref);
           end
           
           
           simTrace = zeros(dim, numSteps + 1);
           simTrace(:,1) = init_state;
           controlTrace = zeros(obj.plant.nI, numSteps + 1);
           controlTrace(:,1) = zeros(obj.plant.nI, 1);
                     
           for i=2:numSteps + 1
               x = simTrace(:, i-1);
               y = obj.plant.C * x;
               y1 = [ref_input; y];
               u = obj.controller.evaluate(y1);
               x0 = [x; u];
               [~,~,x1] = obj.transPlant.initial(x0, obj.controlPeriod);
               n = size(x1, 1);
               x2 = x1(n, :)'; % store computed states to simTrace 
               x2(obj.plant.dim + 1:obj.transPlant.dim) = []; 
               simTrace(:, i) = x2;
               controlTrace(:,i) = u;
           end
           
           obj.simTraces = [obj.simTraces; simTrace];
           obj.controlTraces = [obj.controlTraces; controlTrace];
           simTime = toc(t); 
           
        end
        
        
        % generate a number of simulation traces, used for falsification
        function [sim_Traces, control_Traces, genTime] = generateTraces(obj, init_set, ref_input, numSteps, N)
            % @init_set: initial set of state, a star set
            % @ref_input: reference input, a star set or a vector
            % @numSteps: number of simulation steps
            % @N: number of simulation traces need to be generated
            % @sim_Traces: a cell of simulation traces
            % @control_Traces: a cell of control traces
            % @genTim: generating time
            
            % author: Dung Tran
            % date: 10/2/2019
            
            t = tic;
            
            if ~isa(init_set, 'Star')
                error('Initial set is not a star set');
            end
            
            if init_set.dim ~= obj.plant.dim
                error('Inconsistent dimension between the initial set and the plant');
            end
            
            if ~isa(ref_input, 'Star') && ~isempty(ref_input) && (size(ref_input, 2 ~= 1) && size(ref_input, 1) ~= obj.nI_ref)
                error('Invalid ref_input, should be a star set, empty or a vector with %d elements', obj.nI_ref);
            end
            
            init_states = init_set.sample(N); % sample initial set of states 
                        
            if isempty(ref_input)
                
                M = size(init_states, 2);
                sim_Traces = cell(M, 1); % reset simulation traces
                control_Traces = cell(M, 1); % reset control traces

                for i=1:M
                    [sim_Traces{i}, control_Traces{i}, ~] = obj.simulate(init_states(:,i), ref_input, numSteps);                
                end                
                
            else
                
                if isa(ref_input, 'Star')
                   ref_inputs = ref_input.sample(N); % sample reference inputs
                else % ref_input is a vector
                    ref_inputs = zeros(obj.nI_ref, N);
                    for i=1:N
                        ref_inputs(:, i) = ref_input;
                    end
                end
   
                M = min(size(init_states, 2), size(ref_inputs, 2));
                sim_Traces = cell(M,1); % reset simulation traces
                control_Traces = cell(M,1); % reset control traces
                
                for i=1:M
                    [sim_Traces{i}, control_Traces{i}, ~] = obj.simulate(init_states(:,i), ref_inputs(:,i), numSteps);
                end              
                
            end
            
            obj.simTraces = sim_Traces;
            obj.controlTraces = control_Traces;
            genTime = toc(t);
            

        end
        
        
    end
    
    
    methods % PLOT SIMULATION TRACES
        
        function plotSimTraces(varargin)
            % @index: index of the state needs to be plotted
            % @color: color of trace
            % @marker: marker for plot
            
            % author: Dung Tran
            % date: 10/2/2019
            
            switch nargin
                
                case 2
                    obj = varargin{1};
                    index = varargin{2};
                    color = 'blue';
                    markers = '-x';
                case 3
                    obj = varargin{1};
                    index = varargin{2};
                    color = varargin{3};
                    markers = '-x';
                case 4
                    obj = varargin{1};
                    index = varargin{2};
                    color = varargin{3};
                    markers = varargin{4};
                otherwise
                    error('Invalid number of inputs, should be 1, 2, or 3');
            end
                    
            
            if isempty(obj.simTraces)
                error('simulation traces are empty, please do simulation first, i.e., run simulate method');
            end
            
            n = length(obj.simTraces); % number of simulation traces
            m = length(obj.simTraces{1}); % number of steps
            
            T = 1:m; 
            
            for i=1:n-1
                simTrace = obj.simTraces{i};
                plot(T, simTrace(index, :), markers, 'color', color);
                hold on;                
            end
            
            simTrace = obj.simTraces{n};
            plot(T, simTrace(index, :), markers, 'color', color);
            
        end
        
        
        
    end
    
    
    methods % FALSIFICATION USING SIMULATION
        
        function [safe, falsifyTraces, falsifyTime] = falsify(obj, falsifyPRM)
            % @fasifyPRM: falsification parameters including following
            % inputs:
            %       1) @init_set: initial set of state, a star set
            %       2) @ref_input: reference input, a vector or a star set
            %       3) @numSteps: number of simulation steps 
            %       4) @numTraces: number of traces generated for falsification
            %       5) @controlPeriod: control period
            %       6) @unsafeRegion: a HalfSpace object
            % unsafe region is defined by: U = unsafe_mat * x <= unsafe_vec
            
            % @safe: = 0: unsafe
            %        = 2: unknown (falsification is incomplete)
            % @falsifyTrace: falsified traces
            % @falsifyTime: falsification time
            
            
            % author: Dung Tran
            % date: 10/3/2019
            % update 3/15/2020
            
            t = tic;
            
            initSet = falsifyPRM.init_set;
            ref_input = falsifyPRM.ref_input;
            numSteps = falsifyPRM.numSteps;
            numTraces = falsifyPRM.numTraces;
            obj.controlPeriod = falsifyPRM.controlPeriod;
            unsafe_mat = falsifyPRM.unsafeRegion.G;
            unsafe_vec = falsifyPRM.unsafeRegion.g;
            
            obj.generateTraces(initSet, ref_input, numSteps, numTraces);
                       
            n = length(obj.simTraces);
                        
            m = size(obj.simTraces{1}, 2);
            
            U = HalfSpace(unsafe_mat, unsafe_vec);
            
           obj.falsifyTraces = {}; % a cell of fasified traces
            
            for i=1:n
                simTrace = obj.simTraces{i};
                for j=1:m
                    U.contains(simTrace(:,j));
                    obj.falsifyTraces = [obj.falsifyTraces; simTrace];
                    break;
                end
                
            end
            
            if isempty(obj.falsifyTraces)
                safe = 2;
                fprintf('The safety of the system is unknown, try increase the number of simulation traces to find counter examples');
            else
                safe = 0;
                fprintf('The system is unsafe, %d falsified traces are found', length(obj.falsifyTraces));
            end
            
            
            falsifyTraces = obj.falsifyTraces; 
            falsifyTime = toc(t);
            obj.falsifyTime = falsifyTime;
                       
            
        end
        
        
        
        % PLOT FALSIFICATION TRACES
        function plotFalsifyTraces(varargin)
            % @index: index of the state needs to be plotted
            % @color: color of trace
            % @marker: marker for plot
            
            % author: Dung Tran
            % date: 10/2/2019
            
            switch nargin
                
                case 2
                    obj = varargin{1};
                    index = varargin{2};
                    color = 'blue';
                    markers = '-x';
                case 3
                    obj = varargin{1};
                    index = varargin{2};
                    color = varargin{3};
                    markers = '-x';
                case 4
                    obj = varargin{1};
                    index = varargin{2};
                    color = varargin{3};
                    markers = varargin{4};
                otherwise
                    error('Invalid number of inputs, should be 1, 2, or 3');
            end
                    
            
            if isempty(obj.falsifyTraces)
                error('simulation traces are empty, please do simulation first, i.e., run simulate method or generateTraces method');
            end
            
            n = length(obj.falsifyTraces); % number of falsify traces
            m = size(obj.falsifyTraces{1},2); % number of steps
            
            T = 1:m; 
            
            for i=1:n-1
                falsifyTrace = obj.falsifyTraces{i};
                plot(T, falsifyTrace(index, :), markers, 'color', color);
                hold on;                
            end
            
            falsifyTrace = obj.falsifyTraces{n};
            plot(T, falsifyTrace(index, :), markers, 'color', color);
            
        end
               
        
    end
    
    
    
    
end

