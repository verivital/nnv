classdef new_NNCS < handle
    % Neural network control system class
    % Merging all previous NNCS classes into one
    %   Diego Manzanas - 08/03/2023
    % 
    % todo: add support for neuralODEs (NN with ODEblocklayers) as plant
    %       models for NNCS
    
    properties
        controller = []; % neural-network controller (NN class)
        plant = []; % plant model, could be linear, nonlinear or neural network-based plant (any ODE, HA or neuralODE (NN))
        feedbackMap = []; % a feedback matrix decribes the mapping from a group of outputs of the plant to a group of inputs of the controller
                          
        % neural network control system architecture
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
        execFirst = 'controller'; % default (compute 'plant' or 'controller' reach sets first)
        
        % used for reachable set computation
        ref_I = []; % reference input set
        init_set = []; % initial set for the plant
        reachSetTree = []; % reachable set tree
        reachTime = 0; % reachable set computation time
        controlSet = []; % control signal of the controller over time
        plantSet = []; % reachable set of plant at every control step
        
        numCores = 1; % can do parallel for linear
        controllerReachOptions = struct; % reachMethod = approx-star by default
        plantReachMethod = ''; % by default, direct for linear, zono-lin for nonlinear

        % simulation
        simTraces = {}; % simulation trace
        controlTraces = {}; % control trace

        % falsification
        falsifyTraces = {};
        falsifyTime = 0;

    end
    
    methods % main methods (verify, reachability, ...)
        
        % constructor
        function obj = NNCS(varargin)
            % @controller: a neural net controller
            % @plant: a plant model (a LinearODE, DLinearODE or Neural net)
            % @feedbackMap: a feedback map from outputs of the plant to the
            % input of the controller
            
            switch nargin
                
                case 3
                    controller = varargin{1};
                    plant = varargin{2};
                    feedbackMap = varargin{3};
                case 2
                    controller = varargin{1};
                    plant = varargin{2};
                    feedbackMap = [0];
                otherwise
                    error('Invalid number of arguments');
            end
            
            if ~isa(controller, 'NN') 
                error('The controller is not a valid neural network');
            end
            
            if ~isa(plant, 'LinearODE') && ~isa(plant, 'DLinearODE') && ~isa(plant, 'NonLinearODE') && ~isa(plant, 'DNonLinearODE') && ~isa(plant, 'HybridA')
                % check if it is a neural ode
                % if ~isa(plant, 'NN')
                %     error('The plant is not valid');
                % else
                %     N_l = length(NN.Layers);
                %     ode_status = 0; % check if NN is a neuralode (valid plant, otherwise throw error)
                %     for L=1:N_l
                %         if isa(NN.Layers{L}, 'ODEblockLayer')
                %             ode_status = 1;
                %             break
                %         end
                %     end
                %     if ode_status == 0
                %         error('The plant is not valid');
                %     end
                % end
                error('The plant is not valid');
            end            
                        
            [nO, nI] = size(feedbackMap);
            
            if nI ~= 1
                error('FeedbackMap should have one column');
            end
            
            if isa(plant, 'NN') % if it's not a neuralODE
                obj.nO = plant.OutputSize;
            else
                obj.nO = plant.nO;
            end

            obj.nI = controller.InputSize;

            if nO * obj.nO > obj.nI
                error('Two many feedback inputs');
            end
                        
            obj.controller = controller;
            obj.plant = plant;
            obj.feedbackMap = feedbackMap;
            obj.nI_fb = nO * obj.nO;
            obj.nI_ref = obj.nI - obj.nI_fb;

            % assign default values
            obj.controllerReachOptions.reachMethod = 'approx-star'; % reachMethod = approx-star by default
            % plant
            % if isa(obj.plant, 'NN')
            %     obj.plantReachMethod = struct;
            %     obj.plantReachMethod.reachMethod = 'approx-star';
            % elseif contains(class(obj.plant), 'nonlinear')
            if contains(class(obj.plant), 'nonlinear')
                obj.plantReachMethod = 'lin'; % linearize method for zono or polyzono in CORA
            else
                obj.plantReachMethod = 'direct'; % direct method for linear systems
            end

            
        end
       
        % reachability analysis of nncs using stars (other possible intermediate representations)
        function S = reach(obj, reachPRM)
            % Syntax:
            % S = reach(obj, reachPRM)
            %
            % @reachPRM: reachability parameters including following inputs
            %       1) @init_set: initial set of states of the plant
            %       2) @ref_input: reference input to the controller. it = [] if there is no reference input.
            %       3) @numCores: number of cores used for computation
            %       4) @numSteps: number of reachability steps
            
            % Check input correctness and assign values
            if ~isstruct(reachPRM)
                error("reachPRM must be s struct defining reachability options")
            end
            if ~isfield(reachPRM, "init_set") || ~isfield(reachPRM, "numSteps")
                error("Fields init_set and numSteps are mandatory to perform reachability analysis of nncs.")
            end
            initSet = reachPRM.init_set;
            n_steps = reachPRM.numSteps;
            
            if isfield(reachPRM, "ref_input")
                ref_inputSet = reachPRM.ref_input;
            else
                ref_inputSet = [];
            end
            
            if isfield(reachPRM, 'numCores')
                n_cores = reachPRM.numCores;
            else
                n_cores = 1;
            end
            
            if isfield(reachPRM, 'reachMethod')
                reach_method = reachPRM.reachMethod;
            else
                reach_method = 'approx-star';
            end
            
            % to compute the reachable sets of the plant or controller first
            if isfield(reachPRM, 'execFirst') 
                obj.execFirst = reachPRM.execFirst;
            end
            
            if ~isa(initSet, 'Star')
                error('Initial set of the plant is not a Star');
            end
            
            if ~isempty(ref_inputSet) && ~isa(ref_inputSet, 'Star') && isvector(ref_inputSet)
                ref_inputSet = Star(ref_inputSet, ref_inputSet);
            end
            
            if n_steps < 1
                error('Number of steps should be >= 1');
            end
            
            start_time = tic; 
            obj.reachSetTree = SetTree(n_steps + 1); % initialize reach set tree
            obj.init_set = initSet;
            obj.ref_I = ref_inputSet;
            obj.reachSetTree.addReachSet(initSet, 1); % add the init_set to the reachSetTree
            if contains(class(obj.plant), 'nonlinear')
                obj.plant.cora_set = []; % ensure reach set is empty to start
            end
            obj.plant.intermediate_set = []; % ensure reach set is empty to start
            
            % Create reachOptions for controller
            reachOpt.numCores = n_cores;
            reachOpt.reachMethod = reach_method;
            if n_cores > 1
                reachOpt.reachOption = 'parallel';
            end
            
            for i=2:n_steps + 1
             
                % reachability analysis for  controller
                fb_I = obj.reachSetTree.extract_fb_ReachSet(i - 1);   
                input_set = obj.nextInputSetStar(fb_I{1});
                if i==2 && strcmp(obj.execFirst, 'plant')
                    lb = zeros(obj.plant.nI, 1); % at first step, U = 0 if compute plant first
                    U = Star(lb, lb);
                else
                    U = obj.controller.reach(input_set, reachOpt); % control set at step i   
                end
                if length(U) > 1
                    U1 = Star.get_hypercube_hull(U);
                    U1 = U1.toStar();
                else
                    U1 = U;
                end

                % reachability analysis for plant
                obj.controlSet = [obj.controlSet U1];
                R = obj.plant.stepReachStar(fb_I{1}(length(fb_I{1})), U1);                 
                obj.reachSetTree.addReachSet(R, i);                 
            end
            % Save reach data
            obj.reachTime = toc(start_time);
            S = obj.reachSetTree.flatten();
             
        end

        % verify safety after doing reachability analysis
        % unsafe region defined by: unsafe_mat * x <= unsafe_vec
        function [safe, checkingTime] = check_safety(obj, unsafeRegion, numCores)
            % @unsafeRegion: a Halfpsace object
            %   Unsafe region is defined by: y: unsafe_mat * x <= unsafe_vec
            % @numCores: number of cores using for checking safety
            % @safe: = 1: safe, 0: unsafe or unknown
            
            t = tic;

            unsafe_mat = unsafeRegion.G;
            unsafe_vec = unsafeRegion.g;
            
            [n1, m1] = size(unsafe_mat); 
            [n2, m2] = size(unsafe_vec);
            
            if n1 ~= n2
                error('Inconsistent dimension between unsafe matrix and unsafe vector');
            end
            
            if m1 ~= obj.plant.dim
                error('Inconsistent dimension between unsafe matrix and plant');
            end
            if m2 ~= 1
                error('Invalid unsafe vector');
            end
            
            S = obj.plant.intermediate_reachSet;
            N = length(S);
            j = 0; 
            % set up parallel computing with number of cores (workers)
            if numCores > 1
                obj.start_pool;                
                parfor i=1:N
                    L = S(i).intersectHalfSpace(unsafe_mat, unsafe_vec);                 
                    if ~isempty(L)
                        fprintf('\nThe %d^th reach set reach unsafe region', i);
                        j = j + 1; 
                    end
                end
            else
                for i=1:N
                    L = S(i).intersectHalfSpace(unsafe_mat, unsafe_vec);                 
                    if ~isempty(L)
                        fprintf('\nThe %d^th reach set reach unsafe region', i);
                        j = j + 1; 
                    end
                end
            end
            % return safety result
            if j >= 1
                safe = 0;
            else
                safe = 1;
            end
            
            checkingTime = toc(t);
            
        end
        
        % Verify property (is this _ check_safety necessary?)
        function [safe, counterExamples, verifyTime] = verify(obj, reachPRM, unsafeRegion)
            % Syntax:
            % [safe, counterExamples, verifyTime] = verify(obj, reachPRM, unsafeRegion)
            %
            % @reachPRM: reachability parameters consists of following
            % inputs: 
            %       1) @reachPRM.init_set: initial set
            %       2) @reachPRM.ref_input: reference input
            %       3) @reachPRM.numSteps: number of steps
            %       4) @reachPRM.numCores: number of cores used in reachability analysis
            % @unsafeRegion: a Halfpsace object
            % Usafe region is defined by: y: unsafe_mat * x <= unsafe_vec
            %
            % @safe: = unsafe
            %        = safe
            %        = unknown (due to conservativeness)
            % @counterExamples: an array of star set counterExamples or
            %                   falsified input points
            % @verifyTime: verification time
            
            % author: Dung Tran
            % date:   2/15/2019
            % update: 3/15/2020
            
            t = tic;
            obj.reach(reachPRM);
            falsifyPRM.controlPeriod = obj.plant.controlPeriod;
            falsifyPRM.numSteps = reachPRM.numSteps;
            falsifyPRM.init_set = reachPRM.init_set.getBox;
            if isvector(reachPRM.ref_input)
                falsifyPRM.ref_input = Box(reachPRM.ref_input, reachPRM.ref_input);
            else
                if ~isempty(reachPRM.ref_input)
                    falsifyPRM.ref_input = reachPRM.ref_input.getBox;
                else
                    falsifyPRM.ref_input = [];
                end
            end
            falsifyPRM.unsafeRegion = unsafeRegion;
            falsifyPRM.numTraces = 1000;
            
                        
            [safe1, ~] = obj.check_safety(unsafeRegion.G, unsafeRegion.g, reachPRM.numCores);
            counterExamples = [];
            if safe1 == 1
                safe = 'SAFE';
            else
                [rs, ~, counterExamples, ~, ~, ~] = obj.falsify(falsifyPRM);
                if rs == 1
                    safe = 'UNSAFE';
                else
                    safe = 'UNKNOWN';
                end
            end
            
            verifyTime = toc(t);
            
        end
    
        % live reachability analysis, plot reachable set on the fly and produce video for the analysis
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
    end

    methods % simulation based
        
        % simulate (evaluate) the nncs with specific input and initial state of the plant
        function [simTrace, controlTrace] = evaluate(obj, step, n_steps, x0, ref_input)
            % @step: control step size
            % @N: number of control steps
            % @x0: initial state of the plant
            % @simTrace: simulation trace
            % @controlTrace: control signal correpsonding to simulation trace

            if step <= 0
                error('Invalid control step size');
            end
            
            if n_steps < 1
                error('Invalid number of steps');
            end
            
            if ~isempty(ref_input)
                if size(ref_input, 1) ~= obj.nI_ref
                    error('Inconsistent dimension between reference input vector and number of reference inputs');
                end
            
                if size(ref_input, 2) ~= 1
                    error('Invalid reference input vector');
                end
            end
            
            
%             [~,y1] = obj.plant.evaluate([0 step], x0, 0); % first step simulation
            [~,y1] = obj.plant.evaluate(x0, 0); % first step simulation
            n = size(y1, 1);
            obj.simTraces = [];
            obj.controlTraces = [];
            obj.simTraces = [obj.simTraces y1(n, :)'];
            obj.controlTraces = zeros(obj.controller.nO, 1); % control signal of the first step is zero
      
            if n_steps >= 2
                
                for i=2:n_steps
                    % construct input to the controller
                    l = size(obj.simTraces, 2);
                    m = size(obj.feedbackMap, 1);
                    I = [];
                    for j=1:m
                        if l - obj.feedbackMap(j) <= 0
                            I1 = zeros(obj.plant.nO, 1); 
                            I = [I; I1];
                        else
                            I2 = obj.plant.C * obj.simTraces(:, l - obj.feedbackMap(j));
                            I = [I; I2];
                        end 
                    end
                    
                    I = [ref_input; I];
                    % compute control signal
                    u = obj.controller.evaluate(I);
                    % compute states of the plant                  
%                     [~,y1] = obj.plant.evaluate([0 step], obj.simTrace(:, i-1), u); % first step simulation
                    [~,y1] = obj.plant.evaluate(obj.simTraces(:, i-1), u); % first step simulation
                    n = size(y1, 1);
                    obj.simTraces = [obj.simTraces y1(n, :)']; % store computed states to simTrace                    
                    obj.controlTraces = [obj.controlTraces u]; % store control input to controlTrace
                end
                               
            end
            obj.simTraces = [x0 obj.simTraces]; % add initial state to simtrace            
            simTrace = obj.simTraces;
            controlTrace = obj.controlTraces;
            
        end
        
        % randomly simulate nncs
        function [sim_time, sim_traces, control_traces, sampled_init_states, sampled_ref_inputs] = sample(obj, step, n_steps, init_set, ref_input_set, n_samples)
            % @step: control step
            % @n_steps: number of control steps
            % @init_set: initial state of plant, needed to be a box
            % @ref_input_set: reference input set, needed to be a box
            % @n_samples: number of samples
            % @sim_time: simulation time for n_samples
            % @sim_traces: a cell of simulation traces
            % @control_traces: a cell of control traces

            t = tic; 
            
            if ~isa(init_set, 'Box')
                error('Initial states of the plant should be a box');
            end
            if init_set.dim ~= obj.plant.dim
                error('Inconsistent dimension between initial set of state and plant');
            end
            
            if ~isa(ref_input_set, 'Box')
                error('Reference input set should be a box');
            end
            if ~isempty(ref_input_set) && ref_input_set.dim ~= obj.nI_ref
                error('Inconsitence between reference input set and number of reference inputs in nncs object');
            end
            
            if n_samples < 1
                error('Number of samples shoule be >= 1');
            end
            
            % sampling the network with n_samples of input vector           
            % get sampled input vectors
            X = cell(1, obj.plant.dim);
            V = []; % initial input vectors
            for i=1:obj.plant.dim
                X{1, i} = (init_set.ub(i) - init_set.lb(i)).*rand(n_samples, 1) + init_set.lb(i);
                V = vertcat(V, X{1, i}');
            end
           
            Z = []; % reference input vectors
            
            if ~isempty(ref_input_set)
                Y = cell(1, obj.nI_ref);
                for i=1:obj.nI_ref
                    Y{1, i} = (ref_input_set.ub(i) - ref_input_set.lb(i)).*rand(n_samples, 1) + ref_input_set.lb(i);
                    Z = vertcat(Z, Y{1, i}');
                end                
            end           
            
            sampled_init_states = V;
            sampled_ref_inputs = Z;
            sim_traces = cell(1, n_samples);
            control_traces = cell(1, n_samples);
            
            for i=1:n_samples
                if isempty(Z) % no reference input
                     [sim_traces{1, i}, control_traces{1, i}] = obj.evaluate(step, n_steps, V(:, i), []);
                else
                    [sim_traces{1, i}, control_traces{1, i}] = obj.evaluate(step, n_steps, V(:, i), Z(:, i));
                end
            end
            
            sim_time = toc(t);
            
        end
        
        % automatically falsify nncs using random simulations
        function [falsify_result, falsify_time, counter_sim_traces, counter_control_traces, counter_init_states, counter_ref_inputs] = falsify(obj, step, n_steps, init_set, ref_input_set, unsafeRegion, n_samples)
            % @step: control step size
            % @n_steps: number of control steps
            % @init_set: initial set of the plant, should be a box
            % @ref_input_set: reference input set, should be a box
            % @unsafeRegion: a Halfpsace object
            %   Usafe region is defined by: y: unsafe_mat * x <= unsafe_vec
            % @n_samples: number of simulations used for falsification
            
            % @falsify_result: = 1: counter example exist, = 0: counter
            % example does not exit, -> increase number of samples
            % @falsify_time: falsification time
            % @counter_sim_traces: counter simulation traces
            % @counter_control_traces: counter control traces correpsonding
            % to counter simulation traces
            % @counter_init_states: counter initial states of plant
            % @counter_ref_inputs: counter reference inputs
            
            t = tic; 
            [~, sim_traces, control_traces, sampled_init_states, sampled_ref_inputs] = obj.sample(step, n_steps, init_set, ref_input_set, n_samples);
            
            n = size(sim_traces, 2);
            violate_trace_indexes = [];
            for i=1:n
                violate = NNCS.check_trace(sim_traces{1, i}, unsafeRegion);
                if violate
                    violate_trace_indexes = [violate_trace_indexes i];
                end
            end
            
            if isempty(violate_trace_indexes)
                fprintf('Cannot find counter examples, please consider increasing number of samples for falsification');
                falsify_result = 0;
            else
                fprintf('Counter examples are found, %d traces in %d simluation traces violate safety property', length(violate_trace_indexes), n_samples);
                falsify_result = 1;
            end
            
            n = length(violate_trace_indexes);
            counter_sim_traces = cell(1, n);
            counter_control_traces = cell(1,n);
            counter_init_states = cell(1,n);
            counter_ref_inputs = cell(1,n);
            
            for i=1:n
                counter_sim_traces{1, i} = sim_traces(:, violate_trace_indexes(i));
                counter_control_traces{1, i} = control_traces(:, violate_trace_indexes(i));
                counter_init_states{1, i} = sampled_init_states(:, violate_trace_indexes(i));
                counter_ref_inputs{1, i} = sampled_ref_inputs(:, violate_trace_indexes(i));
            end
            
            falsify_time = toc(t);
           
        end
        
    end
    
    
    methods(Static) % helper methods
        
        % check if a trace violates safety specification
        function violate = check_safety_trace(simTrace, unsafeRegion)
            % @simTrace: a single simulation trace
            % @unsafeRegion: a Halfpsace object
            %   Usafe region is defined by: y: unsafe_mat * x <= unsafe_vec
            % @violate: =1: trace reaches unsafe region
            %           =0: trace does not reach unsafe region
            
            [n, m] = size(simTrace);
            unsafe_mat = unsafeRegion.G;
            unsafe_vec = unsafeRegion.g;
            [n1, m1] = size(unsafe_mat);
            [n2, m2] = size(unsafe_vec);
             
            if n ~= m1
                error('Inconsistent dimension between simTrace and unsafe matrix');
            end
            
            if n1 ~= n2
                error('Inconsistent dimension between unsafe matrix and unsafe vector');
            end
            
            if m2 ~= 1
                error('Invalid unsafe vector, it should have one column');
            end
            
            A = unsafe_mat * simTrace - unsafe_vec; 

            k = 0;
            for i=1:m
                for j=1:n1
                    if A(j, i) <= 0
                        k = k + 1;
                    end
                end
                if k == n1
                    break;
                else
                    k = 0;
                end
            end 
        
            if k == n1
                violate = 1;
            else
                violate = 0;
            end
            
        end
        
    end

    methods % helper methods

        % get next step input set with Stars
        function I = nextInputSetStar(obj, fb_I)
            % @fb_I: feed back input set
            
            l = length(fb_I);
            fb_inputSet = [];
            if l > 0
                for i=1:l
                    if ~isa(fb_I(i), 'Star')
                        error('The %d th feedback input is not a Star', i);
                    end
                    fb_inputSet = [fb_inputSet fb_I(i).affineMap(obj.plant.C, [])];
                end
            end           
                        
            n = size(obj.feedbackMap, 1);          
            
            if isempty(obj.ref_I) && isempty(fb_inputSet)
                I = [];
            end
            if ~isempty(obj.ref_I) && isempty(fb_inputSet)
                new_V = vertcat(obj.ref_I, zeros(obj.nI_fb, obj.ref_I.nVar + 1));
                I = Star(new_V, obj.ref_I.C, obj.ref_I.d);
            end
            
            if ~isempty(obj.ref_I) && ~isempty(fb_inputSet)
                l = length(fb_inputSet);
                nA = fb_inputSet(1).dim;
                I2 = [];
                for i=1:n
                    if l - obj.feedbackMap(i) <= 0
                        I2 = Star(zeros(nA, 1),zeros(nA, 1));
                    else
                        I2 = [I2 fb_inputSet(l - obj.feedbackMap(i))];
                    end                
                end
                I = Star.concatenateStars([obj.ref_I I2]);
            end
            
            if isempty(obj.ref_I) && ~isempty(fb_inputSet)
                l = length(fb_inputSet);
                nA = fb_inputSet(1).dim;
                I2 = [];
                for i=1:n
                    if l - obj.feedbackMap(i) <= 0
                        I2 = Star(zeros(nA, 1),zeros(nA, 1));
                    else
                        I2 = [I2 fb_inputSet(l - obj.feedbackMap(i))];
                    end                
                end
                lb = zeros(obj.nI_ref,1);
                ub = zeros(obj.nI_ref,1);
                I1 = Star(lb,ub);
                I = Star.concatenateStars([I1 I2]);
            end         
            
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
        

%         function S = get_intermediate_reachSet(obj)
%             if isa(obj.plant, 'NonLinearODE')
%                 obj.plant.get_interval_sets();
%                 S = obj.plant.intermediate_reachSet;
%             elseif isa(obj.plant, 'DNonLinearODE')
%                 S = obj.plant.intermediate_reachSet;
%             elseif isa(obj.plant, 'DLinearODE')
%                 S = 
%             elseif isa(obj.plant, 'LinearODE')
%         end
    end

    methods % plot methods 
        % TODO: update, use CORA reach sets and methods to avoid converting all sets to stars if possible (nonlinear)
        
        % plot output reachable set of the nncs
        function plotOutputReachSets(obj, color, map_mat, map_vec)
            % @map_mat: a mapping matrix
            % @map_vec: mapping vector
            % Y = map_mat * X + map_vec
            % @color: color
            
            Y = obj.getOutputReachSet(map_mat, map_vec);
            n = length(Y); % number of control periods
            
            % plot output reach sets
            option = size(map_mat, 1);
            h = obj.plant.reachStep;                       
            if option == 1 % plot 1D, output versus time steps        
                for i=1:n
                    B = Star.get_hypercube_hull(Y{i});
                    ymin = B.lb;
                    ymax = B.ub;
                    y = [ymin ymin ymax ymax ymin];
                    xmin = (i-1)*h;
                    xmax = xmin + h;
                    x = [xmin xmax xmax xmin xmin];
                    plot(x, y, color);
                    hold on;
                end
                T = 0:obj.plant.controlPeriod:n*obj.plant.controlPeriod;
                ax = gca; 
                ax.XTick = T;
            end
            
            if option == 2 || option == 3 % plot 2D or 3D
                G = [];
                for i=1:n
                     G = [G Star.get_hypercube_hull(Y{i})];                   
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
        
        % get output set, specify transformations with map_mat and map_vec
        function Y = getOutputReachSet(obj, map_mat, map_vec)
            % @map_mat: a mapping matrix
            % @map_vec: mapping vector
            % Y = map_mat * X + map_vec
            % @Y: a cell of output reach sets
            
            if isempty(obj.plant.intermediate_reachSet)
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
            n = length(obj.plant.intermediate_reachSet);
            
            Y = cell(1, n);
            for i=1:n
                Y{i} = obj.plant.intermediate_reachSet(i).affineMap(map_mat, map_vec);
            end
            
        end
        
    end
    
end

