classdef NN < handle
    % NN class encodes all type of neural networks supported in NNV
        %   This class is part of the refactoring plan for NNV executed at the
        %   end of 2022. This class generalizes previous constructors for
        %   neuralODEs, FFNNS, CNN, BNN, RNN,...
        %   In addition, we will add a property Connections, in order to support
        %   other type of DL models like residual or U networks.

    % Reachability analysis methods are developed for each individual layer
    
    % Author: Diego Manzanas Lopez
    % Date: 09/28/2022
    % Notes: Code is based on the previous CNN and FFNNS classes written by
    %             Dr. Hoang Dung Tran
    % This is a generalized class, created in the refactoring of NNV in 2022/2023 (NNV 2.0)
    
    properties
        
        Name = 'nn'; % name of the network
        Layers = {}; % An array of Layers, eg, Layers = [L1 L2 ...Ln]
        Connections = []; % A table specifying source and destination layers
%         name2number = containers.Map; % hashmp to match layers name to its number in Layers array (facilitate reach computation graph)
        numLayers = 0; % number of Layers
        numNeurons = 0; % number of Neurons
        InputSize = 0; % number of Inputs
        OutputSize = 0; % number of Outputs
        
        % properties for reach set computation
        reachMethod = 'approx-star';    % reachable set computation scheme, default - 'approx-star'
        relaxFactor = 0; % default - solve 100% LP optimization for finding bounds in 'approx-star' method
        reachOption = []; % parallel option, default - non-parallel computing
        numCores = 1; % number of cores (workers) using in computation
        reachSet = [];  % reachable set for each layers
%         outputSet = [];
        reachTime = []; % computation time for each layers
%         totalReachTime = 0; % total computation time
%         type = 'nn'; % default type is nn (general neural network), other options: 'neuralODE', 'bnn_xnor', 'cnn', 'segmentation', 
        features = {}; % outputs of each layer in an evaluation
        input_vals = {}; % input values to each layer
        input_sets = {}; % input set values for each layer
        dis_opt = []; % display option = 'display' or []
        lp_solver = 'linprog'; % choose linprog as default LP solver for constructing reachable set
        % user can choose 'glpk' or 'linprog' as an LP solver
        
    end
    
    
    methods % constructor, evaluation, sampling, print methods
        
        % constructor
        function obj = NN(varargin)
            
            switch nargin
                % if connections undefined, assume NN is fullyconnected
                % No skipped connections, no multiple connections from any layer
                case 5
                    % parse inputs
                    Layers = varargin{1};
                    cons = varargin{2}; % connections
                    inputsize = varargin{3};
                    outputsize = varargin{4};
                    name = varargin{5};
                    nL = length(Layers); % number of Layers

                    % update object properties
                    obj.Name = name;                % Name of the network
                    obj.Layers = Layers;             % Layers in NN
                    obj.Connections = cons;       % Connections in NN
                    obj.numLayers = nL;             % number of layers
                    obj.InputSize = inputsize;      % input size
                    obj.OutputSize = outputsize; % output size

                case 4
                    % parse inputs
                    Layers = varargin{1};
                    inputsize = varargin{2};
                    outputsize = varargin{3};
                    name = varargin{4};
                    nL = length(Layers); % number of Layers
%                     warning('No connections were specified, we assume each layer is connected to the next layer in the order they were defined in the Layers array');
                    N = length(Layers);
                    sources = 1:N;
                    dests = 2:N+1;
                    conns = table(sources', dests', 'VariableNames', {'Source', 'Destination'});
                    
                    % update object properties
                    obj.Name = name;                % Name of the network
                    obj.Layers = Layers;            % Layers in NN
                    obj.numLayers = nL;             % number of layers
                    obj.InputSize = inputsize;      % input size
                    obj.OutputSize = outputsize;    % output size
                    obj.Connections = conns;

                case 2 % only layers and connections defined
                    Layers = varargin{1};
                    conns = varargin{2};
                    obj.Layers = Layers;
                    obj.Connections = conns;

                case 1 % only layers, assume no sparse connections
                    Layers = varargin{1};
                    N = length(Layers);
                    sources = 1:N;
                    dests = 2:N+1;
                    conns = table(sources', dests', 'VariableNames', {'Source', 'Destination'});
                    obj.Layers = Layers;
                    obj.Connections = conns;

                case 0
                    obj.Layers = {};
                    obj.Connections = [];
                    obj.numLayers = 0;
                    obj.InputSize = 0;
                    obj.OutputSize = 0;

                    
                otherwise
                    error('Invalid number of inputs, should be 0, 1, 2, 3 or 5');
            end
                      
        end
                
        
        % Evaluation of a NN (update this to make use of the connections graph)
        function y = evaluate(obj, x)
            % Evaluation of a NN (compute output of NN given an input)
            % @x: input vector x
            % @y: output vector y
            % @features: output of all layers
            
            for i=1:height(obj.Connections)
                if i > 1
                    x = obj.input_vals{obj.Connections.Source(i)}; % 
                end
                y = obj.Layers{obj.Connections.Source(i)}.evaluate(x);
                % Work on this statement after adding support for residual networks
%                 if isfield(obj.Layers{obj.Connections(i).Destination}, 'input_val')
%                     obj.Layers{obj.Connections(i).Destination}.input_val{end+1} = y; % store layer input (from all sources)
%                 else
%                     obj.Layers{obj.Connections(i).Destination}.input_val = {y}; % store layer input 
%                 end
                % Keep it simple for now 
                if i < height(obj.Connections) % (Last number in destinations is the output)
                    obj.input_vals{obj.Connections.Destination(i)} = y; % store layer input 
                end
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
          
    end
    
    
    methods % reachability analysis method

        function reachOptions = check_reachability_method(obj, reachOptions)
            reach_method = reachOptions.reachMethod;
            if contains(reach_method, "exact")
                for i=1:length(obj.Layers)
                    if isa(obj.Layers{i}, "ODEblockLayer")
                        if ~isa(obj.Layers{i}.odemodel,'LinearODE')
                            warning("Exact reachability is not possible with a neural ODE layer (" + class(obj.Layers{i}.odemodel) + "). Switching to approx-star.");
                            reachOptions.reachMethod = "approx-star";
                        end
                    end
                    if isa(obj.Layers{i}, "SigmoidLayer") || isa(obj.Layers{i}, "TanhLayer")
                        warning("Exact reachability is not possible with layer " + class(obj.Layers{i} + ". Switching to approx-star."));
                        reachOptions.reachMethod = "approx-star";
                    end
                end
            elseif contains(reach_method, "relax-star")
                if ~isfield(reachOptions, "relaxFactor")
                    error("Please, specify a relax factor value to perform relax-star reachability analysis")
                else
                    if reachOptions.relaxFactor < 0 || reachOptions.relaxFactor > 1
                        error('Invalid relaxFactor. The value of relax factor must be between 0 and 1');
                    end
                end
            end
        end
        
        % Update this function to make use of the computational graph
        % (connections) newly introduced. Can access the computed reach set
        % for each layer using the obj.reachSet property (cell array, same
        % size as the layer array

        % Define the reachability function for any general NN
        function outputSet = reach(obj, inputSet, reachOptions)            
            % inputSet: input set (type -> Star, ImageStar, Zono or ImageZono)   
            % reachOptions: reachability options for NN (type -> struct)
            %       Required Fields:
            %           reachMethod (string)
            %               select from {'approx-star','exact-star', 'abs-dom', approx-zono'}
            %       Optional Fields:
            %           numCores (int) -> default = 1
            %           relaxFactor (int) -> default = 0
            %           dis_opt ('display' or []), display options, use for debugging, default = [] (no display)
            %           lp_solver ('glpk', 'linprog') -> default = 'glpk'

            % Parse the inputs

            % TODOs: update display functions

            % Ensure input set is a valid type
            if ~ isa(inputSet,"Star") && ~ isa(inputSet,"ImageStar") && ~ isa(inputSet,"ImageZono") && ~ isa(inputSet,"Zono") 
                error('Wrong input set type. Input set must be of type "Star", "ImageStar", "ImageZono", or "Zono"')
            end

            % Check validity of reachability method
            if exist("reachOptions",'var')
                reachOptions = check_reachability_method(obj, reachOptions);
            else
                reachOptions = struct; % empty options, run with default values
            end

            % Process reachability options
            if ~isstruct(reachOptions)
                error("The reachability parameters must be specified as a struct.")
            else
                if isfield(reachOptions, 'reachMethod')
                    obj.reachMethod = reachOptions.reachMethod;
                end
                if isfield(reachOptions, 'numCores')
                    obj.numCores = reachOptions.numCores;
                end
                if isfield(reachOptions, 'relaxFactor')
                    obj.relaxFactor = reachOptions.relaxFactor;
                end
                if isfield(reachOptions, 'dis_opt')
                    obj.dis_opt = reachOptions.dis_opt; % use for debuging
                end
                if isfield(reachOptions, 'lp_solver')
                    obj.lp_solver = reachOptions.lp_solver;
                end
            end
            
            % Parallel computation or single core?
            if  obj.numCores > 1
                obj.start_pool;
                obj.reachOption = 'parallel';
            end

            % Initialize variables to store reachable sets and computation time
            obj.reachSet = cell(1, height(obj.Connections)); % store input reach sets for each layer
            obj.reachTime = zeros(1, height(obj.Connections)); % store computation time for each connection 
            % Debugging option
            if strcmp(obj.dis_opt, 'display')
                fprintf('\nPerform reachability analysis for the network %s \n', obj.Name);
            end

            % Begin reachability computation
            obj.reachSet{1} = inputSet;
            for i=1:height(obj.Connections)
%                 if strcmp(obj.dis_opt, 'display')
%                     fprintf('\nPerforming analysis for Layer %d (%s)... \n', i-1, obj.Layers{i-1}.Name);
%                 end
                rs = obj.reachSet{obj.Connections.Source(i)}; 
                % Compute reachable set layer by layer
                start_time = tic;
                rs_new = obj.Layers{obj.Connections.Source(i)}.reach(rs, obj.reachMethod, obj.reachOption, obj.relaxFactor, obj.dis_opt, obj.lp_solver);
                obj.reachTime(i) = toc(start_time);
                % Store computed reach set
                if i < height(obj.Connections) % (Last number in destinations is the output)
                    obj.reachSet{obj.Connections.Destination(i)} = rs_new; % store layer input 
                end
%                 obj.reachSet{i-1} = rs_new;
%                 if strcmp(obj.dis_opt, 'display')
%                     fprintf('Reachability analysis for Layer %d (%s) is done in %.5f seconds \n', i-1, obj.Layers{i-1}.Name, obj.reachTime(i-1));
%                     fprintf('The number of reachable sets at Layer %d (%s) is: %d \n', i-1, obj.Layers{i-1}.Name, length(rs_new));
%                 end
            end
            if strcmp(obj.dis_opt, 'display')
                fprintf('Reachability analysis for the network %s is done in %.5f seconds \n', obj.Name, sum(obj.reachTime));
                fprintf('The number ImageStar in the output sets is: %d \n', length(rs_new));
            end
            % Output of the function
%             obj.totalReachTime = sum(obj.reachTime);
            if length(obj.reachSet) < obj.Connections.Destination(end) && length(obj.Layers) >= obj.Connections.Destination(end)
                rs = rs_new; 
                % Compute reachable set layer by layer
                start_time = tic;
                outputSet = obj.Layers{obj.Connections.Destination(end)}.reach(rs, obj.reachMethod, obj.reachOption, obj.relaxFactor, obj.dis_opt, obj.lp_solver);
                obj.reachTime(i) = toc(start_time);
            else
                outputSet = rs_new;
            end

        end
        
    end
    
    methods % clasification and verification methods

        function result = verify_vnnlib(obj, propertyFile, reachOptions)
            % Load specification to verify
            [lb_x, ub_x, property] = load_vnnlib(propertyFile);
            % Create reachability parameters and options
            if contains(reachOptions.reachMethod, "zono")
                X = ImageZono(lb_x', ub_x');
            else
                X = ImageStar(lb_x',ub_x');
            end
            % Compute reachability
            rT = tic;
            Y = obj.reach(X, reachOptions); % Seems to be working
            rT = toc(rT);
            result = verify_specification(Y, property); 
            % Evaluate property
            disp(' ');
            disp('===============================')
            disp('RESULTS')
            disp(' ')
            
            if result == 2
                if contains(net.reachMethod, "exact")
                    result = 0;
                    disp('Property is UNSAT');
                else
                    disp('Property is UNKNOWN');
                end
            elseif result == 1
                disp('Property is SAT');
            else
                disp('Property is UNSAT')
            end
    
            disp("Reachability computation time = "+string(rT) + " seconds")
        end
        
        % Check robustness of output set given a target label
        function [rb, cands] = checkRobust(~, outputSet, correct_id)
            % @outputSet: the outputSet we need to check
            % @correct_id: the correct_id of the classified output
            % @rb: = 1 -> robust
            %      = 0 -> not robust
            %      = 2 -> unknown
            % @cand: possible candidates
            
            R = outputSet;
            
            if contains(class(outputSet),'Image')
                if correct_id > outputSet.numChannel || correct_id < 1
                    error('Invalid correct id');
                end

                R = R.toStar;
            end

            [lb, ub] = R.getRanges;
            [~, max_ub_id] = max(ub);
            cands = [];
            if max_ub_id ~= correct_id
                rb = 2;
                cands = max_ub_id;
            else                   
                
                max_val = lb(correct_id);
                max_cd = find(ub > max_val); % max point candidates
                max_cd(max_cd == correct_id) = []; % delete the max_id

                if isempty(max_cd)
                    rb = 1;
                else            
                    
                    n = length(max_cd);
                    C1 = R.V(max_cd, 2:R.nVar+1) - ones(n,1)*R.V(correct_id, 2:R.nVar+1);
                    d1 = -R.V(max_cd, 1) + ones(n,1)*R.V(correct_id,1);
                    S = Star(R.V, [R.C;C1], [R.d;d1], R.predicate_lb, R.predicate_ub);
                    if S.isEmptySet
                        rb = 2;
                        cands = max_cd;
                    else                       
                        count = 0;
                        for i=1:n
                            if R.is_p1_larger_than_p2(max_cd(i), correct_id)
                                rb = 2;
                                cands = max_cd(i);
                                break;
                            else
                                count = count + 1;
                            end
                        end
                        if count == n
                            rb = 1;
                        end         
                    end    

                end

            end

        end % end check robust
        
        function label_id = classify(obj, inputValue, reachOptions)
            % label_id = classify(inputValue, reachOptions = 'default')
            % inputs
                % @inputValue: input to the NN, typically an image, or an imagestar 
                % reachOptions: 'default' means no reach parameters are
                    % chosen, to be used also when input is not a set
            % @label_id: output index of classified object
            % Assume the network outputs the largest output value

            % parse inputs
            switch nargin
                case 1
                    reachOptions = 'default';
                case 2
                otherwise
                    error('Invalid number of inputs, should be 1 or 2');
            end
            
            % Check if input is not a set, then proceed with classification
            if ~isa(inputValue, 'ImageStar') && ~isa(inputValue, 'Star') && ~isa(inputValue, 'Zono') && ~isa(inputValue, 'ImageZono')
                y = obj.evaluate(in_image);
                y = reshape(y, [obj.OutputSize, 1]);
                [~, label_id] = max(y); 
            else
                % For input sets, compute reach set and then classify
                if ~isstruct(reachOptions)
                    reachOptions.method = 'approx-star'; % default reachability method
                    reachOptions.numOfCores = 1;
                end
                RS = obj.reach(inputValue, reachOptions);
                % For now, we are converting the set to ImageStar and then
                % computing max values for classification
                n = length(RS);
                label_id = cell(n, 1);
                for i=1:n
                    rs = RS(i);
                    if isa(rs, 'Star') || isa(rs, 'Zono')
                        rs = rs.toImageStar(rs.dim, 1, 1);
                    elseif isa(rs, 'ImageZono')
                        rs = rs.toImageStar;
                    end
                    new_rs  = ImageStar.reshape(rs, [obj.OutputSize(1) 1 1]);                    
                    max_id = new_rs.get_localMax_index([1 1], [obj.OutputSize(1) 1], 1);
                    label_id{i} = max_id(:, 1);
                end
            end
        end       
    end
    
end

