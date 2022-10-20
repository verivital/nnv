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
    %             Dung Tran
    % This is a generalized class, created after the refactoring of NNV in 2022
    
    properties
        
        Name = 'nn'; % name of the network
        Layers = {}; % An array of Layers, eg, Layers = [L1 L2 ...Ln]
        Connections = []; % A table specifying source and destination layers
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
%         reachTime = []; % computation time for each layers
%         totalReachTime = 0; % total computation time
%         type = 'nn'; % default type is nn (general neural network), other options: 'neuralODE', 'bnn_xnor', 'cnn', 'segmentation', 
        features = {}; % outputs of each layer in an evaluation
        dis_opt = []; % display option = 'display' or []
        lp_solver = 'glpk'; % choose linprog as default LP solver for constructing reachable set
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
                    warning('No connections were specified, we assume each layer is connected to the next layer in the order they were defined in the Layers array');

                    % update object properties
                    obj.Name = name;                % Name of the network
                    obj.Layers = Layers;            % Layers in NN
                    obj.numLayers = nL;             % number of layers
                    obj.InputSize = inputsize;      % input size
                    obj.OutputSize = outputsize;    % output size

                case 2 % only layers and connections defined
                    Layers = varargin{1};
                    conns = varargin{2};
                    obj.Layers = Layers;
                    obj.Connections = conns;

                case 0
                    obj.Layers = {};
                    obj.Connections = [];
                    obj.numLayers = 0;
                    obj.InputSize = 0;
                    obj.OutputSize = 0;

                    
                otherwise
                    error('Invalid number of inputs, should be 0, 4 or 5');
            end
                      
        end
                
        
        % Evaluation of a NN
        function y = evaluate(obj, x)
            % Evaluation of a NN (compute output of NN given an input)
            % @x: input vector x
            % @y: output vector y
            % @features: output of all layers
            
            y = x;
            for i=1:obj.numLayers
                y = obj.Layers{i}.evaluate(y); 
                obj.features{i} = y;
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

            % Ensure input set is a valid type
            if ~ isa(inputSet,"Star") && ~ isa(inputSet,"ImageStar") && ~ isa(inputSet,"ImageZono") && ~ isa(inputSet,"Zono") 
                error('Wrong input set type. Input set must be of type "Star", "ImageStar", "ImageZono", or "Zono"')
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
            obj.reachSet = cell(1, obj.numLayers);
            obj.reachTime = zeros(1, obj.numLayers);
            % Debugging option
            if strcmp(obj.dis_opt, 'display')
                fprintf('\nPerform reachability analysis for the network %s...', obj.Name);
            end
            % Begin reachability computation
            rs = inputSet;
            for i=2:obj.numLayers+1
                if strcmp(obj.dis_opt, 'display')
                    fprintf('\nPerforming analysis for Layer %d (%s)...', i-1, obj.Layers{i-1}.Name);
                end
                % Compute reachable set layer by layer
                start_time = tic;
                rs_new = obj.Layers{i-1}.reach(rs, obj.reachMethod, obj.reachOption, obj.relaxFactor, obj.dis_opt, obj.lp_solver);
                obj.reachTime(i-1) = toc(start_time);
                % Update input set to the next layer
                rs = rs_new;
                % Store computed reach set
                obj.reachSet{i-1} = rs_new;
                if strcmp(obj.dis_opt, 'display')
                    fprintf('\nReachability analysis for Layer %d (%s) is done in %.5f seconds', i-1, obj.Layers{i-1}.Name, obj.reachTime(i-1));
                    fprintf('\nThe number of reachable sets at Layer %d (%s) is: %d', i-1, obj.Layers{i-1}.Name, length(rs_new));
                end
            end
            if strcmp(obj.dis_opt, 'display')
                fprintf('\nReachability analysis for the network %s is done in %.5f seconds', obj.Name, sum(obj.reachTime));
                fprintf('\nThe number ImageStar in the output sets is: %d', length(rs_new));
            end
            % Output of the function
            obj.totalReachTime = sum(obj.reachTime);
            outputSet = rs_new;

        end
        
    end
    
    methods % clasification and verification methods
        
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
        
        % To generalize things like this, it may be a little more
        % complicated, have to double check if we could do counterexamples
        % with any sets and any NN type
        %%%%%%%%%%%%%%%%%%%%%%%
%  THIS NEEDS FIXING, HOW TO DO COUNTER EXAMPLES WITH ANY SET AND ANY NN 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [robust, counterExamples] = verifyRobustness(obj, inputValue, targetID, reachOptions)
            % [robust, counterExamples] = verifyRobustness(obj, inputValue, targetID, reachOptions)
            % @robust: = 1: the network is robust
            %          = 0: the network is notrobust
            %          = 2: robustness is uncertain
            % @counterExamples: a set of counter examples 
            % @inputValue = input to NN, i.e. ImageStar, Star...
            % @targetID = label (class) target of the inputValue
            % @reachOptions = reachability options for the robustness analysis
            
            switch nargin
                case 3
                    reachOptions.method = 'approx-star';
                    reachOptions.numOfCores = 1;
                case 4
                    if ~isstruct(reachOptions)
                        error('Reachability options must be defined as a struct.')
                    end
                otherwise
                    error('Invalid number of inputs, should be 2, 3, or 4');
            end
                        
            label_id = obj.classify(inputValue, reachOptions);      
            n = length(label_id); 
            % check the correctness of classifed label
            counterExamples = [];
            incorrect_id_list = []; 
            for i=1:n
                ids = label_id{i};
                m = length(ids);
                id1 = [];
                for j=1:m
                    if ids(j) ~= targetID
                       id1 = [id1 ids(j)];
                    end
                end
                incorrect_id_list = [incorrect_id_list id1];             
                
                % construct counter example set
                if ~isempty(id1) % Check for incorrect classifications
                    % if it's not a set, return the example input
                    if ~isa(inputValue, 'ImageStar') && ~isa(inputValue, 'Star') && ~isa(inputValue, 'Zono') && ~isa(inputValue, 'ImageZono')
                        counterExamples = in_image; 
                    elseif strcmp(reachOptions.method, 'exact-star') && isa(inputValue, 'ImageStar') % Can only return counter examples if method used was exact
                        rs = obj.outputSet(i);
                        L = length(id1); 
                        for l=1:L                    
                            [new_C, new_d] = ImageStar.addConstraint(rs, [1, 1, correct_id], [1, 1, id1(l)]);  
                            counter_IS = ImageStar(in_image.V, new_C, new_d, in_image.pred_lb, in_image.pred_ub);
                            counterExamples = [counterExamples counter_IS];
                        end
                    end
                end
            end
            
            
            if isempty(incorrect_id_list)
                robust = 1;
                fprintf('\n=============================================');
                fprintf('\nTHE NETWORK IS ROBUST');
                fprintf('\nClassified index: %d', correct_id);
                
            else
                
                if strcmp(method, 'exact-star')
                    robust = 0;
                    fprintf('\n=============================================');
                    fprintf('\nTHE NETWORK IS NOT ROBUST');
                    fprintf('\nLabel index: %d', correct_id);
                    fprintf('\nClassified index:');
                    
                    n = length(incorrect_id_list);
                    for i=1:n
                        fprintf('%d ', incorrect_id_list(i));
                    end
                    
                else
                    robust = 2;
                    fprintf('\n=============================================');
                    fprintf('\nThe robustness of the network is UNCERTAIN due to the conservativeness of approximate analysis');
                    fprintf('\nLabel index: %d', correct_id);
                    fprintf('\nPossible classified index: ');
                                       
                    n = length(incorrect_id_list);
                    for i=1:n
                        fprintf('%d ', incorrect_id_list(i));
                    end                   
                    fprintf('\nPlease try to verify the robustness with exact-star (exact analysis) option');
                end
                
            end
            
        end
        
    end
    
end

