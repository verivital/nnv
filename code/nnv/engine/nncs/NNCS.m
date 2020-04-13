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
        controlSet = []; % control signal of the controller over time
        simTrace = []; % simulation trace
        controlTrace = []; % control trace
    end
    
    methods
        
        %constructor
        function obj = NNCS(varargin)
            % @controller: a neural net controller
            % @plant: a plant model (a LinearODE, DLinearODE or Neural net)
            % @feedbackMap: a feedback map from outputs of the plant to the
            % input of the controller
            
            % author: Dung Tran
            % date: 11/1/2018
            % update: 3/14/2020
            
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
            
            
            if ~isa(controller, 'FFNN') && ~isa(controller, 'FFNNS')
                error('The controller is not a feedforward neural network');
            end
            
            if ~isa(plant, 'LinearODE') && ~isa(plant, 'DLinearODE') && ~isa(plant, 'NonLinearODE') && ~isa(plant, 'DNonLinearODE')
                error('The plant is not a linear or nonlinear ode system in discrete or continuous time');
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
                       
        % reachability analysis of nncs using stars
        % output reach set of controller is a single star
        % the plant reachable set is a zonotope
        function [P, reachTime] = reach(obj, reachPRM)
             % @reachPRM: reachability parameters including following inputs
             %       1) @init_set: initial set of states of the plant
             %       2) @ref_input: reference input to the controller. it = [] if there is no reference input.
             %       4) @numCores: number of cores used for computation
             %       5) @numSteps: number of reachability steps

             % author: Dung Tran
             % date: 11/16/2018
             
             start_time = tic; 
             initSet = reachPRM.init_set;
             ref_inputSet = reachPRM.ref_input;
             n_cores = reachPRM.numCores;
             n_steps = reachPRM.numSteps;
             
             if ~isa(initSet, 'Star')
                 error('Initial set of the plant is not a Star');
             end

             if ~isempty(ref_inputSet) && ~isa(ref_inputSet, 'Star')
                 error('The reference input set is not a Star');
             end

             if n_steps < 1
                 error('Number of steps should be >= 1');
             end
             
             obj.reachSetTree = SetTree(n_steps + 1); % initialize reach set tree
             obj.init_set = initSet;
             obj.ref_I = ref_inputSet;
             obj.reachSetTree.addReachSet(initSet, 1); % add the init_set to the reachSetTree
             
             for i=2:n_steps + 1
                 
                 % reachability analysis for  controller
                 fprintf('Reachability analysis for the controller \n');
                 fb_I = obj.reachSetTree.extract_fb_ReachSet(i - 1);   
                 input_set = obj.nextInputSetStar(fb_I{1});
                 [U,~] = obj.controller.reach(input_set, 'exact-star', n_cores); % control set at step i
                 U1 = Star.get_hypercube_hull(U);   
                 
                 % reachability analysis for plant
                 fprintf('\nReachability analysis for the plant \n');
                 U1 = U1.toStar();
                 obj.controlSet = [obj.controlSet U1];
                 R = obj.plant.stepReachStar(fb_I{1}(length(fb_I{1})), U1);                 
                 obj.reachSetTree.addReachSet(R, i);                 
             end
             reachTime = toc(start_time);
             obj.reachTime = reachTime;
             P = obj.reachSetTree.flatten();
             obj.totalNumOfReachSet = obj.reachSetTree.getTotalNumOfReachSet();
             
             
        end

        % get next step input set with Stars
        function I = nextInputSetStar(obj, fb_I)
            % @fb_I: feed back input set
            
            % author: Dung Tran
            % date: 11/18/2018
            
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
                        
                        P = Polyhedron('lb', zeros(nA, 1), 'ub', zeros(nA, 1));
                        
                        I2 = [I2 Conversion.toStar(P)];
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
                        P = Polyhedron('lb', zeros(nA, 1), 'ub', zeros(nA, 1));
                        I2 = [I2 Conversion.toStar(P)];
                    else

                        I2 = [I2 fb_inputSet(l - obj.feedbackMap(i))];

                    end                

                end
                
                lb = zeros(obj.nI_ref,1);
                ub = zeros(obj.nI_ref,1);
                I1 = Polyhedron('lb', lb, 'ub', ub);
                I = Star.concatenateStars([Conversion.toStar(I1) I2]);
                
            end         
            
        end
        
        
        % verify safety after doing reachability analysis
        % unsafe region defined by: unsafe_mat * x <= unsafe_vec
        function [safe, checkingTime] = check_safety(obj, unsafe_mat, unsafe_vec, numOfCores)
            % @unsafe_mat: unsafe region matrix
            % @unsafe_vec: unsafe region vector
            % @numOfCores: number of cores using for checking safety
            % @safe: = 1: safe, 0: unsafe or unknown
            
            % author: Dung Tran
            % date: 1/18/2019
            
            t = tic;
            
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
            if numOfCores > 1
                poolobj = gcp('nocreate'); % If no pool, do not create new one.
                if isempty(poolobj)
                    parpool('local', numOfCores); 
                else
                    if poolobj.NumWorkers ~= numOfCores
                        delete(poolobj); % delete the old poolobj
                        parpool('local', numOfCores); % start the new one with new number of cores
                    end                    
                end                   
                             
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
            
            if j >= 1
                safe = 0;
            else
                safe = 1;
            end
            
            checkingTime = toc(t);
            
        end
        
        % simulate (evaluate) the nncs with specific input and initial state of the plant
        function [simTrace, controlTrace] = evaluate(obj, step, n_steps, x0, ref_input)
            % @step: control step size
            % @N: number of control steps
            % @x0: initial state of the plant
            % @simTrace: simulation trace
            % @controlTrace: control signal correpsonding to simulation
            % trace
            
            % author: Dung Tran
            % date: 1/29/2019
            
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
            
            
            [~,y1] = obj.plant.evaluate([0 step], x0, 0); % first step simulation
            n = size(y1, 1);
            obj.simTrace = [];
            obj.controlTrace = [];
            obj.simTrace = [obj.simTrace y1(n, :)'];
            obj.controlTrace = zeros(obj.controller.nO, 1); % control signal of the first step is zero
      
            if n_steps >= 2
                
                for i=2:n_steps
                    
                    % construct input to the controller
                    l = size(obj.simTrace, 2);
                    m = size(obj.feedbackMap, 1);
                    I = [];
              
                    for j=1:m
              
                        if l - obj.feedbackMap(j) <= 0
                            I1 = zeros(obj.plant.nO, 1); 
                            I = [I; I1];
                        else
                            I2 = obj.plant.C * obj.simTrace(:, l - obj.feedbackMap(j));
                            I = [I; I2];

                        end 

                    end
                    
                    I = [ref_input; I];
                   
                    % compute control signal
                    u = obj.controller.evaluate(I);
                    % compute states of the plant                  
                    [~,y1] = obj.plant.evaluate([0 step], obj.simTrace(:, i-1), u); % first step simulation
                    n = size(y1, 1);
                    obj.simTrace = [obj.simTrace y1(n, :)']; % store computed states to simTrace                    
                    obj.controlTrace = [obj.controlTrace u]; % store control input to controlTrace
                end
                               
            end
            obj.simTrace = [x0 obj.simTrace]; % add initial state to simtrace            
            simTrace = obj.simTrace;
            controlTrace = obj.controlTrace;
            
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
            
            % author: Dung Tran
            % date: 1/31/2019
            
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
        function [falsify_result, falsify_time, counter_sim_traces, counter_control_traces, counter_init_states, counter_ref_inputs] = falsify(obj, step, n_steps, init_set, ref_input_set, unsafe_mat, unsafe_vec, n_samples)
            % @step: control step size
            % @n_steps: number of control steps
            % @init_set: initial set of the plant, should be a box
            % @ref_input_set: reference input set, should be a box
            % @unsafe_mat: unsafe matrix
            % @unsafe_vec: unsafe vector
            % @n_samples: number of simulations used for falsification
            
            % @falsify_result: = 1: counter example exist, = 0: counter
            % example does not exit, -> increase number of samples
            % @falsify_time: falsification time
            % @counter_sim_traces: counter simulation traces
            % @counter_control_traces: counter control traces correpsonding
            % to counter simulation traces
            % @counter_init_states: counter initial states of plant
            % @counter_ref_inputs: counter reference inputs
                        
            
            % author: Dung Tran
            % date: 1/31/2019
            
            t = tic; 
            [~, sim_traces, control_traces, sampled_init_states, sampled_ref_inputs] = obj.sample(step, n_steps, init_set, ref_input_set, n_samples);
            
            n = size(sim_traces, 2);
            violate_trace_indexes = [];
            for i=1:n
                violate = NNCS.check_trace(sim_traces{1, i}, unsafe_mat, unsafe_vec);
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
    
    
    methods(Static)
        
        % parse a SpaceEx file for reachability analysis
        function plant = parseSpaceEx(xmlFile,rootID,outputName)
            % Parse a file from SpaceEx and load a plant into NNV
            % 
            % plant = parseSpaceEx(xmlFile,rootID,outputName)
            % 
            % Inputs:
            %   xmlFile - path to the xml-file that contains the SpaceEx-model
            %   rootID -  ID of SpaceEx component to be used as root component
            %           (specified as a string)
            %   outputName - name of the generated CORA model (specified as string)
            %            (CORA model will be saved in current working directory)
            % 
            % plant - NNV plant model obtained from the CORA file
           
            % -------- Example ---------
            % plant = NNCS.parsePlant('platoon_continuous.xml','platoon11', 'platoonC');
            % 
            % See also: spaceex2cora.m
            spaceex2cora(xmlFile, rootID, outputName, pwd);
            fi = char(outputName);
            fi = [fi '.m'];
            fid = fopen(fi);
            if fid == -1, error('Cannot open file'); end
            filetext = fileread(fi);
            matches = regexp(filetext,'dynamics =','match');
            if length(matches) > 1, error('We do not support HA modeling as part of the NNCS class'); end
            while true
                tline = fgetl(fid);
                if length(tline) > 9
                    if string(tline(1:8)) == "dynamics"
                        dynamics = 'nonlinear';
                        disp('Creating the nonlinear plant');
                        break
                    end
                end
                if length(tline) > 5
                   if string(tline(1:4)) == "dynA"
                       dynamics = 'linear';
                       disp('Creating the linear plant');
                       break
                   end
                end
            end
            if dynamics == "nonlinear"
                fclose(fid);
                str = split(string(tline),' = ');
                str = str(2);
                str = split(str,'(');
                dynamics = str(1);
                str = str(2);
                str = split(str,',');
                num_states = str2double(str(1));
                num_inputs = str2double(str(2));
                str = split(str(3),'@');
                strfunc = str2func(str(2));
                if contains(dynamics,'DT')
                    Ts = 0.1;
                    C = eye(num_states);
                    plant = DNonLinearODE(num_states, num_inputs,strfunc,Ts,C);
                    disp('Plant created with default values.');
                    disp(plant);
                    disp('Please see DNonLinearODE for more details.');
                else
                    reachStep = 0.01;
                    controlPeriod = 0.1;
                    C= eye(num_states);
                    plant = NonLinearODE(num_states,num_inputs,strfunc,reachStep,controlPeriod,C);
                    disp('Plant created with default values.');
                    disp(plant);
                    disp('Please see NonLinearODE for more details.');
                end
            elseif dynamics == "linear"
                % Initialize matrices for linear model
                dynA = [];
                dynB = [];
                dynC = [];
                dynD = [];
                dynRow = 1;
                dynCol = 1;
                p = 1;
                % Obtain matrix A
                while true
                    tline = fgetl(fid);
                    if string(tline) == "-1"
                        error('Stop')
                    end
                    if  contains(string(tline),"dynB")
                        disp('Moving to matrix B');
                        break
                    end
                    str = split(string(tline),'[');
                    str2 = str(end);
                    str = split(str2,';');
                    if str2 == str
                        tmpl = split(str,',');
                        for p = dynCol:dynCol+length(tmpl)-1
                            if contains(string(tmpl(p+1-dynCol)),"]")
                                tt = split(string(tmpl(p+1-dynCol)),']');
                                if ~isnan(str2double(tt(1)))
                                    dynA(dynRow,p) = str2double(tt(1));
                                end
                                break
                            elseif contains(string(tmpl(p+1-dynCol)),"...")
                                tt = split(string(tmpl(p+1-dynCol)),'...');
                                if ~isnan(str2double(tt(1)))
                                    dynA(dynRow,p) = str2double(tt(1));
                                end
                            elseif ~isnan(str2double(tmpl(p+1-dynCol)))
                                dynA(dynRow,p) = str2double(tmpl(p+1-dynCol));
                            end
                        end
                        dynCol = p;
                    else
                        for i=dynRow:dynRow+length(str)-1
                            tmpl = split(str(i+1-dynRow),',');
                            for k=dynCol:dynCol+length(tmpl)-1
                                if contains(string(tmpl(k+1-dynCol)),"]")
                                    tt = split(string(tmpl(k+1-dynCol)),']');
                                    if ~isnan(str2double(tt(1)))
                                        dynA(i,p) = str2double(tt(1));
                                    end
                                    break                                    
                                elseif contains(string(tmpl(k+1-dynCol)),"...")
                                    tt = split(string(tmpl(k+1-dynCol)),'...');
                                    if ~isnan(str2double(tt(1)))
                                        dynA(dynRow,p) = str2double(tt(1));
                                    end
                                elseif ~isnan(str2double(tmpl(k+1-dynCol)))
                                    dynA(i,k) = str2double(tmpl(k+1-dynCol));
                                end
                            end
                            dynCol = 1;
                        end
                        dynRow = i;
                    end
                end
                % Obtain matrix B
                dynRow = 1;
                dynCol = 1;
                p = 1;
                while true
                    tline = fgetl(fid);
                    if string(tline) == "-1"
                        error('Stop')
                    end
                    if  contains(string(tline),"dync") || contains(string(tline),"dynC")
                        disp('Moving to matrix C');
                        break
                    end
                    str = split(string(tline),'[');
                    str2 = str(end);
                    str = split(str2,';');
                    if str2 == str
                        tmpl = split(str,',');
                        for p = dynCol:dynCol+length(tmpl)-1
                            if contains(string(tmpl(p+1-dynCol)),"]")
                                tt = split(string(tmpl(p+1-dynCol)),']');
                                if ~isnan(str2double(tt(1)))
                                    dynB(dynRow,p) = str2double(tt(1));
                                end
                                break
                            elseif contains(string(tmpl(p+1-dynCol)),"...")
                                tt = split(string(tmpl(p+1-dynCol)),'...');
                                if ~isnan(str2double(tt(1)))
                                    dynB(dynRow,p) = str2double(tt(1));
                                end
                            elseif ~isnan(str2double(tmpl(p+1-dynCol)))
                                dynB(dynRow,p) = str2double(tmpl(p+1-dynCol));
                            end
                        end
                        dynCol = p;
                    else
                        for i=dynRow:dynRow+length(str)-1
                            tmpl = split(str(i+1-dynRow),',');
                            for k=dynCol:dynCol+length(tmpl)-1
                                if contains(string(tmpl(k+1-dynCol)),"]")
                                    tt = split(string(tmpl(k+1-dynCol)),']');
                                    if ~isnan(str2double(tt(1)))
                                        dynB(i,p) = str2double(tt(1));
                                    end
                                    break                                    
                                elseif contains(string(tmpl(k+1-dynCol)),"...")
                                    tt = split(string(tmpl(k+1-dynCol)),'...');
                                    if ~isnan(str2double(tt(1)))
                                        dynB(dynRow,p) = str2double(tt(1));
                                    end
                                elseif ~isnan(str2double(tmpl(k+1-dynCol)))
                                    dynB(i,k) = str2double(tmpl(k+1-dynCol));
                                end
                            end
                            dynCol = 1;
                        end
                        dynRow = i;
                    end
                end
                dynRow = 1;
                dynCol = 1;
                p = 1;
                % Obtaining matrix C
                while true
                    tline = fgetl(fid);
                    if string(tline) == "-1"
                        error('Stop')
                    end
                    if  contains(string(tline),"dynamics")
                        disp('Finishing the matrices');
                        break
                    end
                    str = split(string(tline),'[');
                    str2 = str(end);
                    str = split(str2,';');
                    if str2 == str
                        tmpl = split(str,',');
                        for p = dynCol:dynCol+length(tmpl)-1
                            if contains(string(tmpl(p+1-dynCol)),"]")
                                tt = split(string(tmpl(p+1-dynCol)),']');
                                if ~isnan(str2double(tt(1)))
                                    dynC(dynRow,p) = str2double(tt(1));
                                end
                                break
                            elseif contains(string(tmpl(p+1-dynCol)),"...")
                                tt = split(string(tmpl(p+1-dynCol)),'...');
                                if ~isnan(str2double(tt(1)))
                                    dynC(dynRow,p) = str2double(tt(1));
                                end
                            elseif ~isnan(str2double(tmpl(p+1-dynCol)))
                                dynC(dynRow,p) = str2double(tmpl(p+1-dynCol));
                            end
                        end
                        dynCol = p;
                    else
                        for i=dynRow:dynRow+length(str)-1
                            tmpl = split(str(i+1-dynRow),',');
                            for k=dynCol:dynCol+length(tmpl)-1
                                if contains(string(tmpl(k+1-dynCol)),"]")
                                    tt = split(string(tmpl(k+1-dynCol)),']');
                                    if ~isnan(str2double(tt(1)))
                                        dynC(i,p) = str2double(tt(1));
                                    end
                                    break                                    
                                elseif contains(string(tmpl(k+1-dynCol)),"...")
                                    tt = split(string(tmpl(k+1-dynCol)),'...');
                                    if ~isnan(str2double(tt(1)))
                                        dynC(dynRow,p) = str2double(tt(1));
                                    end
                                elseif ~isnan(str2double(tmpl(k+1-dynCol)))
                                    dynC(i,k) = str2double(tmpl(k+1-dynCol));
                                end
                            end
                            dynCol = 1;
                        end
                        dynRow = i;
                    end
                end
                try
                    dynD = zeros(size(dynC,1),size(dynB,2));
                    plant = LinearODE(dynA,dynB,dynC,dynD);
                catch
                    try
                        dynC = dynC';
                        dynD = zeros(size(dynC,1),size(dynB,2));
                        plant = LinearODE(dynA,dynB,dynC,dynD);
                    catch
                        assignin('base','dynA',dynA);
                        assignin('base','dynB',dynB);
                        assignin('base','dynC',dynC');
                        assignin('base','dynD',dynD');
                        disp('All matrices have been saved to the workspace as dynA, dynB, dynC and dynD');
                        error('We have found inconsistency in the matrix dimensions.')
                    end
                end
                disp(plant);
                fclose(fid);
            end
        end
        
        % check if a trace violates safety specification
        function violate = check_trace(simTrace, unsafe_mat, unsafe_vec)
            % @simTrace: a single simulation trace
            % @unsafe_mat: unsafe matrix to specify unsafe region
            % @unsafe_vec: unsafe vector to specify unsafe region:
            % unsafe_mat * x <= unsafe_vec
            % @violate: =1: trace reaches unsafe region
            %           =0: trace does not reach unsafe region
            
            [n, m] = size(simTrace);
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
    
    methods
        % verify method
        function [safe, counterExamples, verifyTime] = verify(obj, reachPRM, unsafeRegion)
            % @reachPRM: reachability parameters consists of following
            % inputs: 
            %       1) @reachPRM.init_set: initial set
            %       2) @reachPRM.ref_input: reference input
            %       3) @reachPRM.numSteps: number of steps
            %       4) @reachPRM.numCores: number of cores used in reachability analysis
            %       5) @reachPRM.plantTimStep: reachability step for the plant
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
            % date:   2/15/2019
            % update: 3/15/2020
            
            
            
            
        end
    end
    
end

