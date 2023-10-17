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
    % This is a generalized class, created in the refactoring of NNV in 2023 (NNV 2.0)
    %    It supports FFNNS, CNN, SEGNET, BNN and RNN from previous NNV version

    % Update: Neelanjana Pal
    % Date: 07/07/2023
    
    properties
        
        Name = 'nn'; % name of the network
        Layers = {}; % An array of Layers, eg, Layers = [L1 L2 ...Ln]
        Connections = []; % A table specifying source and destination layers
        numLayers = 0; % number of Layers
        numNeurons = 0; % number of Neurons
        InputSize = 0; % number of Inputs
        OutputSize = 0; % number of Outputs
        
        % properties for reach set  and evaluation computation
        reachMethod = 'approx-star';    % reachable set computation scheme, default - 'approx-star'
        relaxFactor = 0; % default - solve 100% LP optimization for finding bounds in 'approx-star' method
        reachOption = []; % parallel option, default - non-parallel computing
        numCores = 1; % number of cores (workers) using in computation
        reachSet = {};  % reachable set for each layers
        reachTime = []; % computation time for each layers
        features = {}; % outputs of each layer in an evaluation
        input_vals = {}; % input values to each layer (cell array of cells of layer input values)
        input_sets = {}; % input set values for each layer (cell array of cells of layer input sets)
        dis_opt = []; % display option = 'display' or []
        lp_solver = 'linprog'; % choose linprog as default LP solver for constructing reachable set user can choose 'glpk' or 'linprog' as an LP solver
        
        % To facilitate graph computation flow
        name2indx = []; % Match name to index in nnvLayers list
    end


    methods % main methods (constructor, evaluation, reach)
        
        % constructor
        function obj = NN(varargin)
            
            switch nargin
                % if connections undefined, assume NN is fullyconnected
                % No skipped connections, no multiple connections from any layer
                case 5
                    % parse inputs
                    Layers = varargin{1};
                    conns = varargin{2}; % connections
                    inputSize = varargin{3};
                    outputSize = varargin{4};
                    name = varargin{5};
                    nL = length(Layers); % number of Layers

                case 4
                    % parse inputs 
                    Layers = varargin{1};
                    if isa(varargin{2}, 'table')
                        conns = varargin{2};
                        inputSize = varargin{3};
                        outputSize = varargin{4};
                        name = 'nn';
                    else
                        conns = [];
                        inputSize = varargin{2};
                        outputSize = varargin{3};
                        name = varargin{4};
                    end
                    nL = length(Layers); % number of Layers

                case 2 % only layers and connections defined
                    Layers = varargin{1};
                    conns = varargin{2};
                    nL = length(Layers);
                    name = 'nn';
                    inputSize = 0;
                    outputSize = 0;

                case 1 % only layers, assume each layer is only connected to the next one
                    Layers = varargin{1};
                    nL = length(Layers);
                    conns = [];
                    name = 'nn';
                    inputSize = 0;
                    outputSize = 0;

                case 0
                    name = 'nn';
                    Layers = {};
                    conns = [];
                    nL = 0;
                    inputSize = 0;
                    outputSize = 0;

                otherwise
                    error('Invalid number of inputs, should be 0, 1, 2, 3 or 5');
            end

            % update object properties
            obj.Name = name;              % Name of the network
            obj.Layers = Layers;          % Layers in NN
            obj.Connections = conns;       % Connections in NN
            obj.numLayers = nL;           % number of layers
            obj.InputSize = inputSize;    % input size
            obj.OutputSize = outputSize;  % output size
                      
        end
        
        % Evaluation of a NN (test it, fix for neuralode with multiple outputs)
        function y = evaluate(obj, x)
            % Evaluate NN given an input sample
            % y = NN.evaluate(x)
            % @x: input vector x
            % @y: output vector y
            
            % Two options to exectute evaluation
            % 1) Connections are defined
            if ~isempty(obj.Connections)
                 y = obj.evaluate_withConns(x);
            % 2) No connections defined, execute each layer consecutively
            else
                y = obj.evaluate_noConns(x);
            end

        end

        function y = evaluateSequence(obj, x)
            % Evaluation of a NN (compute output of NN given an input)
            % @x: input vector x
            % @y: output vector y
            
            for i=1:height(obj.Connections)
                if i > 1
                    x = obj.input_vals{1,i}; % 
                end

                L = obj.Layers{i};
                y = L.evaluateSequence(x);
                % Keep it simple for now 
                if i < height(obj.Connections) % (Last number in destinations is the output)
                    obj.input_vals{1,i+1} = y; % store layer input 
                end
            end
        end
        
        % evaluate parallel
        function y = evaluate_parallel(obj, inputs)
            % @inputs: array of inputs
            % @y: output vector
            
            n = length(inputs);
            y = cell(1, n);
            
            if obj.numCores < 1
                error("Invalid number of Cores");
            elseif obj.numCores == 1
                for i=1:n
                    y{i} = obj.evaluate(inputs{i});
                end
            else
                obj.start_pool();
                parfor i=1:n
                    y{i} = obj.evaluate(inputs{i});
                end
            end
        end

        % Define the reachability function for any general NN (test it)
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
            if ~isa(inputSet,"Star") && ~isa(inputSet,"ImageStar") && ~isa(inputSet, "VolumeStar")...
                    && ~isa(inputSet,"ImageZono") && ~isa(inputSet,"Zono") 
                error('Wrong input set type. Input set must be of type "Star", "ImageStar", "VolumeStar", "ImageZono", or "Zono"')
            end

            % Check validity of reachability method
            if exist("reachOptions",'var')
                reachOptions = validate_reach_options(obj, reachOptions);
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
                else
                    obj.numCores = 1;
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

            % Debugging option
            if strcmp(obj.dis_opt, 'display')
                fprintf('\nPerform reachability analysis for the network %s \n', obj.Name);
            end
            
            % Perform reachability based on connections or assume no skip/sparse connections
            if isempty(obj.Connections)
                outputSet = obj.reach_noConns(inputSet);
            else 
                outputSet = obj.reach_withConns(inputSet);
            
            end

        end
    
        function outputSet = reachSequence(obj, inputSet, reachOptions)
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

            % Check validity of reachability method
            if exist("reachOptions",'var')
                reachOptions = validate_reach_options(obj, reachOptions);
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

            % Debugging option
            if strcmp(obj.dis_opt, 'display')
                fprintf('\nPerform reachability analysis for the network %s \n', obj.Name);
            end
            
            outputSet = obj.reach_withConns(inputSet);
            
        end
    
    end
    
    
    methods % secondary methods (verification, safety, robustness...)
        
        % Verify a VNN-LIB specification
        function result = verify_vnnlib(obj, propertyFile, reachOptions)
            
            % Load specification to verify
            property = load_vnnlib(propertyFile);
            lb = property.lb;
            ub = property.ub;

            % Create reachability parameters and options
            if contains(reachOptions.reachMethod, "zono")
                X = ImageZono(lb, ub);
            else
                X = ImageStar(lb,ub);
            end

            % Compute reachability
            Y = obj.reach(X, reachOptions); % Seems to be working
            result = verify_specification(Y, property.prop); 

            % Modify in case of unknown and exact
            if result == 2
                if contains(obj.reachMethod, "exact")
                    result = 0;
                end
            end
    
        end
        
        % Check robustness of output set given a target label
        function rb = checkRobust(obj, outputSet, target)
            % rb = checkRobust(~, outputSet, target)
            % 
            % @outputSet: the outputSet we need to check
            %
            % @target: the correct_id of the classified output or a halfspace defining the robustness specification
            % @rb: = 1 -> robust
            %      = 0 -> not robust
            %      = 2 -> unknown
            
            % Process set
            if ~isa(outputSet, "Star")
                nr = length(outputSet);
                R = Star;
                for s=1:nr
                    R = outputSet(s).toStar;
                end
            else
                R = outputSet;
            end
            % Process robustness spec
            if ~isa(target, "HalfSpace")
                if isscalar(target)
                    target = obj.robustness_set(target, 'max');
                else
                    error("Target must be a HalfSpace (unsafe/not robust region) or as a scalar determining the output label target.")
                end
            end
            % Check robustness
            rb = verify_specification(R, target);

        end % end check robust
        
        % Verify robutness of a NN given an input set and target (id or HalfSpace)
        function result = verify_robustness(obj, inputSet, reachOptions, target)
            % Compute reachable set
            R = obj.reach(inputSet, reachOptions);
            % Check robustness
            result = obj.checkRobust(R, target);
            % Modify in case of unknown and exact
            if result == 2
                if contains(obj.reachMethod, "exact")
                    result = 0;
                end
            end
        end
        
        % Classify input set (one or more label_id?)
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
        
        % Get robustness bound of a NN wrt to an input and disturbance value
        function robustness_bound = get_robustness_bound(varargin)
            % Find maximum robustness value, i.e., maximum disturbance bound
            % that the network is still robust to (iterative search)
            %
            % ---- Syntax ----
            %     [robustness_bound, t] = NN.get_robustness_bound(input_vec, init_dis_bound, tol, max_steps, lb_allowable, ub_allowable, un_robust_reg, reachOptions, n_samples);
            %
            % ---- Inputs ----
            % 1: @input_vec: input point
            % 2: @init_dis_bound: initial disturbance bound
            % 3: @dis_bound_step: a step to increase/decrease disturbance bound
            % 4: @max_steps: maximum number of steps for searching
            % 5: @lb_allowable: allowable lower bound of disturbed input:    lb_allowable(i) <= x'[i]
            % 6: @ub_allowable: allowable upper bound of disturbed output:    ub_allowable(i) >= x'[i] 
            %         x' is the disturbed vector by disturbance bound, |x' - x| <= dis_bound
            %         x'[i] >= lb_allowable, x'[i] <= ub_allowable[i]
            % 7: un_robust_region:  (halfspace(s))
            % 7: @reachOptions (see NN.reach)
            % 8: @n_samples: number of samples used to find counter examples
            % 
            % ---- Outputs ----
            % @robustness_bound: robustness bound w.r.t @input_vec
            %
            
            % author: Diego Manzanas Lopez
            % date: 02/14/2023
            
            switch nargin
                
                case 9
                    obj = varargin{1}; % NN object
                    input_vec = varargin{2}; % input vec
                    init_dis_bound = varargin{3}; % initial disturbance bound for searching
                    tolerance = varargin{4}; % tolerance (accuracy) for searching
                    max_steps = varargin{5}; % maximum searching steps
                    lb_allowable = varargin{6}; % allowable lower bound on disturbed inputs
                    ub_allowable = varargin{7}; % allowable upper bound on disturbed inputs
                    un_robust_reg = varargin{8}; % un-robust region 
                    reachOptions = varargin{9}; % reachability analysis options
                    
                case 7
                    obj = varargin{1}; % FFNNS object
                    input_vec = varargin{2}; % input vec
                    init_dis_bound = varargin{3}; % initial disturbance bound for searching
                    tolerance = varargin{4}; % tolerance (accuracy) for searching
                    max_steps = varargin{5}; % maximum searching steps
                    lb_allowable = []; % allowable lower bound on disturbed inputs
                    ub_allowable = []; % allowable upper bound on disturbed inputs
                    un_robust_reg = varargin{6}; % un-robust region 
                    reachOptions = varargin{7}; % reachability analysis method 
                    
                otherwise
                    error('Invalid number of input arguments (should be 6 or 8).');
            end
            
            % Check reach options
            if ~isstruct(reachOptions)
                error("Please define reachOptions. See NN.reach for indications.")
            elseif ~isfield(reachOptions, 'reachMethod')
                obj.reachMethod = 'approx-star'; % default
            else
                obj.reachMethod = reachOptions.reachMethod;
            end

            % Search for disturbance bound
            k = 1;
            b = init_dis_bound;
            bmax = 0;
            while (k < max_steps)
                % Create input set
                inputSet = obj.create_input_set(input_vec, b, lb_allowable, ub_allowable);
                % Reach set
                outputSet = obj.reach(inputSet, reachOptions);
                % Check robustness
                robust = obj.checkRobust(outputSet, un_robust_reg);

                if robust == 1
                    bmax = b;
                    b = b + tolerance;
                else
                    b = b - tolerance;
                    if b == bmax || b <= 0
                        break; % finish function here
                    end
                end
                k = k + 1;
            end

            % Return result
            if k == max_steps
                if bmax == 0
                    robustness_bound = [];
                else
                    disp("Bound search stopped at kmax, but it is not guaranteed that max bound is found. Try incrementing " + ...
                        "the number of steps to search for maximal distrubance bound");
                    robustness_bound = bmax;
                end
            else
                if bmax == 0
                    robustness_bound = [];
                else
                    robustness_bound = bmax;
                end
            end
            
        end
        
        % Check safety of NN and generate counter examples
        function [result, counter_inputs] = verify_safety(obj, I, U, reachOptions, n_samples)
            %
            % ---- Syntax ----
            % [result, counter_inputs] = NN.verify_safety(I, U, reachOptions, n_samples)
            %
            % ---- Inputs ----
            % 1: @I: input set, need to be a star set
            % 2: @U: unsafe region, a set of HalfSpaces
            % 3: @reachOptions: = 'star' -> compute reach set using stars
            % 4: @n_samples : number of simulations used for falsification if using over-approximate reachability analysis
            % note: n_samples = 0 -> do not do falsification
            %
            % ---- Outputs ----
            % @result: = 1 -> safe, = 0 -> unsafe, = 2 -> uncertain 
            % @counter_inputs
            
            % author: Diego Manzanas 
            % date: 02/10/2023
            
            % Ensure unsafe region is defined            
            if isempty(U)
                error('Please specify unsafe region using Half-space class');
            end
            
            % performing reachability analysis
            R = obj.reach(I, reachOptions);        
                      
            % censure correct set format for checking safety
            n = length(R);
            m = length(U);
            if ~isa(R, 'Star')
                R1 = Star;
                for i=1:n
                    R1(i) = R(i).toStar; % transform to star sets
                end
            else
                R1 = R;
            end
            
            % Possible counter inputs
            violate_inputs = [];

            % Check for possible safety violation
            if isfield(reachOptions, 'numCores')
                if reachOptions.numCores > 1
                    parfor i=1:n
                        for j=1:m
                            S = R1(i).intersectHalfSpace(U(j).G, U(j).g);
                            if ~isempty(S) && strcmp(method, 'exact-star')
                                I1 = Star(I.V, S.C, S.d); % violate input set
                                violate_inputs = [violate_inputs I1];
                            else
                                violate_inputs = [violate_inputs S];
                            end
    
                        end
                    end
                end
            else
                for i=1:n
                    for j=1:m
                        S = R1(i).intersectHalfSpace(U(j).G, U(j).g);
                        if ~isempty(S) && strcmp(reachOptions.reachMethod, 'exact-star')
                            I1 = Star(I.V, S.C, S.d); % violate input set
                            violate_inputs = [violate_inputs I1];
                        else
                            violate_inputs = [violate_inputs S];
                        end
                    end
                end
            end
            
            % Generate counter examples of approx method, otherwise return unsafe input regions
            if isempty(violate_inputs)
                result = 1;  
                counter_inputs = []; 
            else
                if strcmp(reachOptions.reachMethod, 'exact-star')
                    result = 0;  
                    counter_inputs = violate_inputs; % exact-method return complete counter input set
                else
                    if n_samples == 0
                        result = 2;
                        counter_inputs = [];
                    else
                        counter_inputs = obj.falsify(I, U, n_samples);
                        if isempty(counter_inputs)
                            result = 2;
                        else
                            result = 0;
                        end
                    end
                end
            end
        end
        
        % Generate falsification traces given a property, NN and input set
        function counter_inputs = falsify(obj, I, U, n_samples)
            %
            % ---- Syntax ----
            % counter_inputs = NN.falsify(I, U, n_samples)
            %
            % ---- Inputs ----
            % I = input set
            % U: unsafe region (HalfSpace)
            % n_samples: numer of samples used for falsification
            % 

            % init output var
            counter_inputs = [];

            % check U is properly defined
            if ~isa(U, "HalfSpace")
                error('Unsafe region (U) is not a HalfSpace.')
            end

            % Check input set
            if ~isa(I, "ImageZono") && ~isa(I, "ImageStar") && ~isa(I, "Star") && ~isa(I, "Zono") && ~isa(I, "Box")
                error("Input set (I) must be a ImageZono, ImageStar, Star, Zono or Box.")
            end

            % Transform to Star to generate samples (could also implement sample methods there in the future, it may be faster)
            if isa(I, "ImageZono")
                for i=1:length(I)
                    I(i) = I(i).toImageStar;
                end
            elseif isa(I, "Zono")
                for i=1:length(I)
                    I(i) = I(i).toStar;
                end
            end
            
            % Begin falsification atempts
            % Check n_samples
             if n_samples < 1
                error('Invalid number of samples. n_samples > 0');
            end
            % Generate sample inputs
            V = I.sample(n_samples);
            n = size(V, 2); % number of samples 
            m = length(U); % number of HalfSpaces (unsafe/not robust regions)
            % Evaluate and add to solution if it reaches the unsafe region U
            for i=1:n
                y = obj.evaluate(V(:, i));
                for j=1:m
                    if U(j).contains(y)
                        counter_inputs = [counter_inputs V(:, i)];
                        break
                    end
                end
            end
        end

    end % end verification methods


    methods % specific NN methods (sequence, classification, semantic seg...)
        
        % verify robustness of semantic segmentation tasks (CAV2021)
        function [riou, rv, rs, n_rb, n_mis, n_unk, n_att, ver_rs, eval_seg_ims] = verify_segmentation(obj, in_images, ground_truths, reachOptions)
            % 
            % --- Syntax ----
            % [riou, rv, rs, n_rb, n_mis, n_unk, n_att, ver_rs, eval_seg_ims] = verify(obj, in_images, ground_truths, reachOptions)
            %
            % ---- Inputs ----
            % @in_images: an array of input set(s)
            % @ground_truths: an array of ground truth images (input images without attack)
            % @reachOptions: reachability parameters
            %
            % ---- Outputs ----
            % @riou: robust iou
            % @rv: robustness value (percentage of correctly classified pixels)
            % @rs: robustness sensitivity (ratio of (misclassified + unknown pixels)/attacked pixels)
            % @n_rb: number of robust pixels
            % @n_mis: number of misclassified pixels
            % @n_unk: number of unknown pixels
            % @n_att: number of attacked pixels
            % @ver_rs: verified output reachable set, used for plot
            % @eval_seg_ims: evaluated output labels, used for plot
            
            % Process inputs
            n1 = length(ground_truths);
            n2 = length(in_images);
            if n1 ~= n2 || n1 ~= 1
                error("Inputs must be one image as input (ImageStar) and a 1x1 cell containing ground truth labels");
            end
            
            % Initialize output variables
            riou = zeros(1,n1);
            rv = zeros(1, n1);
            rs = zeros(1, n1);
            n_rb = zeros(1, n1);
            n_mis = zeros(1, n1);
            n_unk = zeros(1, n1);
            n_att = zeros(1, n1);
            ver_rs = cell(1, n1);
                       
            % compute reachable set
            obj.reach(in_images, reachOptions);
            Rout = obj.reachSet{end-1}; % reachset prior to output layer
            
            % compute ground truth output segmentation image
            eval_seg_ims = obj.evaluate_parallel(ground_truths);
            
            % compute number of correctly classified pixels
            if obj.OutputSize == 0 % unspecified
                outS = obj.Layers{end}.OutputSize;
                n_pixels = outS(1) * outS(2);
            else % specfied when defining NN
                n_pixels = obj.OutputSize(1) * obj.OutputSize(2);
            end
                        
            % obtain verify rearch set
            % 'unknown': the label of the pixel is unknown, may be correct may be not
            % 'misclass': the label of the pixel is misclassified 
            
            % we introduce 2 more classes for the reachable set
            % 1) unknown class
            % 2) misclassification class
            numClasses = length(obj.Layers{end}.Classes);
                      
            if obj.numCores > 1
                parfor i=1:n1

                    gr_seg_im = eval_seg_ims{i};
                    pc_rs = obj.getPixelClassReachSet(Rout, i);
                    [h, w] = size(pc_rs);
                    ver_im = zeros(h, w);
                    n_mis_ct = 0;
                    n_unk_ct = 0;                  
                    n_att(i) = in_images(i).getNumAttackedPixels; 

                    for j=1:h
                        for k=1:w
                            pc = pc_rs{j, k};
                            if size(pc, 1) == 1
                                if pc == gr_seg_im(j,k)
                                    ver_im(j,k) = pc; %(robust pixel)
                                else
                                    ver_im(j,k) = numClasses; % misclass (unrobust pixel)
                                    n_mis_ct = n_mis_ct + 1;
                                end
                            else
                                c = (pc == gr_seg_im(j, k));
                                if sum(c) ~= 0
                                    n_unk_ct = n_unk_ct + 1;
                                    ver_im(j,k) = numClasses - 1; % unknown pixel
                                else
                                    ver_im(j,k) = numClasses; % unrobust pixel
                                    n_mis_ct = n_mis_ct + 1;
                                end
                                
                            end
      
                        end
                    end                
                    ver_rs{i} = ver_im;
                    n_mis(i) = n_mis_ct;
                    n_unk(i) = n_unk_ct;
                    n_rb_ct = n_pixels - n_mis_ct - n_unk_ct;
                    n_rb(i) = n_rb_ct;
                    iou = jaccard(gr_seg_im, ver_im);
                    iou = iou(~isnan(iou));
                    riou(i) = sum(iou)/length(iou);
                    rv(i) = n_rb_ct/n_pixels;
                    rs(i) = (n_mis_ct + n_unk_ct)/n_att(i);
                end
            else
                 for i=1:n1

                    gr_seg_im = eval_seg_ims{i};
                    pc_rs = obj.getPixelClassReachSet(Rout, i);
                    [h, w] = size(pc_rs);
                    ver_im = zeros(h, w);
                    n_mis_ct = 0;
                    n_unk_ct = 0;                  
                    n_att(i) = in_images(i).getNumAttackedPixels;

                    for j=1:h
                        for k=1:w
                            pc = pc_rs{j, k};
                            if size(pc, 1) == 1
                                if pc == gr_seg_im(j,k)
                                    ver_im(j,k) = pc; %(robust pixel)
                                else
                                    ver_im(j,k) = numClasses; % misclass (unrobust pixel)
                                    n_mis_ct = n_mis_ct + 1;
                                end
                            else
                                c = (pc == gr_seg_im(j, k));
                                if sum(c) ~= 0
                                    n_unk_ct = n_unk_ct + 1;
                                    ver_im(j,k) = numClasses - 1; % unknown pixel
                                else
                                    ver_im(j,k) = numClasses; % unrobust pixel
                                    n_mis_ct = n_mis_ct + 1;
                                end
                                
                            end
      
                        end
                    end                
                    ver_rs{i} = ver_im;
                    n_mis(i) = n_mis_ct;
                    n_unk(i) = n_unk_ct;
                    n_rb_ct = n_pixels - n_mis_ct - n_unk_ct;
                    n_rb(i) = n_rb_ct;
                    iou = jaccard(gr_seg_im, ver_im);
                    iou = iou(~isnan(iou));
                    riou(i) = sum(iou)/length(iou);
                    rv(i) = n_rb_ct/n_pixels;
                    rs(i) = (n_mis_ct + n_unk_ct)/n_att(i);
                end
            end
            
        end
        
        % verify robustness of classification RNNs (HSCC2023)
        function result = verify_sequence_robustness(obj, input_seq, epsilon, target, reachOptions)
            %
            % ---- Syntax ----
            % result = NN.verify_sequence_robustness(inputs, epsilon, target, reachOptions)
            %
            % ---- Inputs ----
            % inputs:
            %
            % ---- Outputs ---- 
            % @result: = 1: the network is robust
            %          = 0: the network is notrobust
            %          = 2: robustness is uncertain
            
            % Process inputs
            x = input_seq;
            n = size(x, 2);
            eps = epsilon;
            y  = target; % label/idx target
            ny = size(y,2);
            
            % Prepare input for reachability analysis
            X = [];
            for i=1:n
                X = [X Star(x(:,i) - eps, x(:, i) + eps)]; % construct sequence of input sets
            end
                        
            % compute reach sets
            Y = obj.reach(X, reachOptions);
            
            % Verification results
            if ny == n
                result = zeros(n,1);
                for i=1:n
                    result(i) = obj.checkRobust(Y(i), target(i));
                end
            else
                result = obj.checkRobust(Y(end), target);
            end
            
        end 

    end % end methods
    

    % helper functions
    methods
        
        % Check reachability options defined are allowed
        function reachOptions = validate_reach_options(obj, reachOptions)
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
                        warning("Exact reachability is not possible with layer " + class(obj.Layers{i}) + ". Switching to approx-star.");
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
            % Ensure reachability is done in a single core
            if ~contains(reach_method, 'exact')
                reachOptions.reachOption = 'single';
                reachOptions.numCores = 1;
            end
        end

        % Create input set based on input vector and bounds
        function R = create_input_set(obj, x_in, disturbance, lb_allowable, ub_allowable) % assume tol is applied to every vale of the input
            % R = create_input_set(obj, x_in, disturbance, lb_allowable, ub_allowable)
            % @R = Star or ImageStar

            lb = x_in;
            ub = x_in;
            n = numel(x_in); % number of elements in array
            % Apply disturbance
            lb = lb - disturbance;
            ub = ub + disturbance;
            % Check is disturbance value is allowed (lb)
            if ~isempty(lb_allowable)
                for i=1:n
                    if lb(i) < lb_allowable(i)
                        lb(i) = lb_allowable(i);
                    end
                end
            end
            % Check is disturbance value is allowed (ub)
            if ~isempty(ub_allowable)
                for i=1:n
                    if ub(i) > ub_allowable(i)
                        ub(i) = ub_allowable(i);
                    end
                end
            end
            
            R = obj.set_or_imageset(lb, ub);

        end
        
        % Given input bounds, create input set
        function R = set_or_imageset(obj, lb, ub)
            if length(ub) ~= numel(ub) % not a vector, so ImageStar
                if contains(obj.reachMethod, 'zono')
                    R = ImageZono(lb,ub);
                else
                    R = ImageStar(lb,ub);
                end
            else
                if contains(obj.reachMethod, 'zono')
                    R = Zono(lb,ub);
                else
                    R = Star(lb,ub);
                end
                for k=1:length(obj.Layers)
                    if isa(obj.Layers{k}, 'ImageInputLayer') || isa(obj.Layers{k}, "Conv2DLayer") || contains(class(obj.Layers{k}), "Pooling")
                        if contains(obj.reachMethod, 'zono')
                            R = ImageZono(lb,ub);
                        else
                            R = ImageStar(lb,ub);
                        end
                        break;
                    end
                end
            end
        end

        % Create unsafe/not robust region from a target label of a classification NN
        function Hs = robustness_set(obj, target, class_type)
            % @Hs: unsafe/not robust region defined as a HalfSpace
            %  - target: label idx of the given input set
            %  - class_type: assume max, but could also be min like in ACAS Xu ('min', 'max')

            outSize = obj.OutputSize;
            if outSize == 0 % output size was not properly defined when creating the NN
                layer = obj.Layers{end};
                if isa(layer, "FullyConnectedLayer")
                    outSize = length(layer.Bias);
                else
                    error("Output size is set to 0, but it must be >= 1");
                end
            elseif target > outSize
                error("Target idx must be less than or equal to the output size of the NN.");
            end

            % Define HalfSpace Matrix and vector
            G = ones(outSize,1);
            G = diag(G);
            G(target, :) = [];
            if strcmp(class_type, "max")
                G = -G;
                G(:, target) = 1;
            elseif strcmp(class_type, "min")
                G(:, target) = -1;
            end
%             g = zeros(height(G),1);

            % Create HalfSapce to define robustness specification
            Hs = [];
            for i=1:height(G)
                Hs = [Hs; HalfSpace(G(i,:), 0)];
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
        
        % evaluate NN when no connections are defined
        function y = evaluate_noConns(obj, x)
            y = x;
            % reset eval related parameters
            obj.features = cell(1, length(obj.Layers));
            for i=1:obj.numLayers
                y = obj.Layers{i}.evaluate(y);
                obj.features{i} = y;
            end
        end
        
        % evaluate NN based on connections table (test it)
        function y = evaluate_withConns(obj, x)
            % rest eval values
            obj.features = cell(1, obj.numLayers); % store output for each layer
            obj.input_vals = cell(1, obj.numLayers); % store input for each layer
           
            % Evaluate layer-by-layer based on connection graph
            for i=1:height(obj.Connections)
                % 1) Get name and index of layer
                source = obj.Connections.Source(i);
                % ensure we get just the name and not specific properties
                source = split(source, '/');
                source = source{1};
                source_indx = obj.name2indx(source); % indx in Layers array
                
                % 2) Get input to layer
                if i > 1
                    x = obj.input_vals{source_indx}; % get inputs to layer unless it's the first layer
                else
                    obj.input_vals{1} = x;
                end
                
                % 3) evaluate layer
                exec_len = length(obj.features);
                % ensure layer has not been evaluated yet
                if isempty(obj.features{source_indx})
                    if isa(obj.Layers{source_indx}, 'MaxUnpooling2DLayer')
                        obj.Layers{source_indx}.MaxPoolIndx = obj.Layers{obj.name2indx(obj.Layers{source_indx}.PairedMaxPoolingName)}.MaxIndx;
                        obj.Layers{source_indx}.MaxPoolSize = obj.Layers{obj.name2indx(obj.Layers{source_indx}.PairedMaxPoolingName)}.InputSize;
                    end
                    y = obj.Layers{source_indx}.evaluate(x);
                    obj.features{source_indx} = y;
                else

                end
                
                % 4) save inputs to destination layer
                dest = obj.Connections.Destination(i);
                dest = split(dest, '/');
                dest_name = dest{1};
                dest_indx = obj.name2indx(dest_name); % indx in Layers array
                % check if there source has multiple inputs (concat layer, unpooling ...)
                if length(dest) > 2
                    % unpooling option
                    if isa(obj.Layers{dest_indx}, 'MaxUnpooling2DLayer')
                        destP = dest{2};
                        if strcmp(destP, 'in')
                            obj.input_vals{dest_indx} = y; % store layer input
                        else
                            error("Destination not valid ("+string(obj.Connections.Destination(i))+")");
                        end
                    % concatenation layers (2 or more inputs with specific order {in1, in2, ... inX})
                    elseif contains(class(obj.Layers{dest_indx}), 'Concatenation')
                        destP = dest{2};
                        if startsWith(destP, 'in')
                            input_num = str2double(destP(3:end));
                        else
                            error("Destination not valid ("+string(obj.Connections.Destination(i))+")");
                        end
                        obj.input_vals{dest_indx}{input_num} = y;
                    % concatenation layers (2 or more inputs with specific order {in1, in2, ... inX})
                    elseif isa(obj.Layers{dest_indx}, 'DepthConcatenationLayer')
                        destP = dest{2};
                        if startsWith(destP, 'in')
                            input_num = str2double(destP(3:end));
                        else
                            error("Destination not valid ("+string(obj.Connections.Destination(i))+")");
                        end
                        obj.input_vals{dest_indx}{input_num} = y;
                    elseif isa(obj.Layers{dest_indx}, 'AdditionLayer')
                        destP = dest{2};
                        if startsWith(destP, 'in')
                            input_num = str2double(destP(3:end));
                        else
                            error("Destination not valid ("+string(obj.Connections.Destination(i))+")");
                        end
                        obj.input_vals{dest_indx}{input_num} = y;
                    else
                        error("Destination not valid ("+string(obj.Connections.Destination(i))+")");
                    end
                else
                    obj.input_vals{dest_indx} = y;
                end
            end
                
            % Check if last layer is executed (default is last layer is not executed, but in the 
            % case of the PixelClassificationLayer this is necessary)
            % Assume last layer in array is the output layer
            if isempty(obj.features{end})
                x = obj.input_vals{end};
                y = obj.Layers{end}.evaluate(x);
                obj.features{end} = y;
            end
        end
        
        % reach NN when no connections are defined (test it)
        function outSet = reach_noConns(obj, inSet)
            % Initialize variables
            obj.reachSet = cell(1, obj.numLayers);
            obj.reachTime = zeros(1, obj.numLayers);
            if strcmp(obj.dis_opt, 'display')
                fprintf('\nPerform reachability analysis for the network %s...', obj.Name);
            end
            % Begin reach set computation
            rs = inSet;
            for i=2:obj.numLayers+1
                if strcmp(obj.dis_opt, 'display')
                    fprintf('\nPerforming analysis for Layer %d (%s)...', i-1, obj.Layers{i-1}.Name);
                end
                start_time = tic;
                rs_new = obj.Layers{i-1}.reach(rs, obj.reachMethod, obj.reachOption, obj.relaxFactor, obj.dis_opt, obj.lp_solver);
                obj.reachTime(i-1) = toc(start_time); % track reach time for each layer
                rs = rs_new; % get input set to next layer
                obj.reachSet{i-1} = rs_new; % store output set for layer
            end
            % Output
            outSet = rs_new;
        end

        % reach NN based on connections table (test it)
        function outSet = reach_withConns(obj, inSet)
            % Initialize variables to store reachable sets and computation time
            obj.reachSet = cell(1, obj.numLayers);
            obj.reachTime = zeros(1, obj.numLayers);
            obj.input_sets = cell(1, height(obj.Connections)); % store input reach sets for each layer
            if strcmp(obj.dis_opt, 'display')
                fprintf('\nPerform reachability analysis for the network %s...', obj.Name);
            end

            % Begin reachability computation
            for i=1:height(obj.Connections)

                % 1) Get name and index of layer
                source = obj.Connections.Source(i);
                % ensure we get just the name and not specific properties
                source = split(source, '/');
                source = source{1};
%                 fprintf('\n layer reachability: %s', source);
                source_indx = obj.name2indx(source); % indx in Layers array
                
                % 2) Get input to layer
                if i > 1
                    inSet = obj.input_sets{source_indx}; % get inputs to layer unless it's the first layer
                end
                
                % 3) reach layer
                exec_len = length(obj.reachSet);
                % ensure layer has not been eexcuted yet
                if exec_len >= source_indx && isempty(obj.reachSet{source_indx})
                    if strcmp(obj.dis_opt, 'display')
                        fprintf('\nPerforming analysis for Layer %d (%s)...', i-1, source);
                    end
                    t = tic;
                    if isequal(class(obj.Layers{1,1}), 'SequenceInputLayer')
                        outSet = obj.Layers{source_indx}.reachSequence(inSet, obj.reachMethod, obj.reachOption, obj.relaxFactor, obj.dis_opt, obj.lp_solver);
                    else
                        outSet = obj.Layers{source_indx}.reach(inSet, obj.reachMethod, obj.reachOption, obj.relaxFactor, obj.dis_opt, obj.lp_solver);
                    end
                    obj.reachTime(source_indx) = toc(t);
                    obj.reachSet{source_indx} = outSet;
                end
                
                % 4) save inputs to destination layer
                dest = obj.Connections.Destination(i);
                dest = split(dest, '/');
                dest_name = dest{1};
                dest_indx = obj.name2indx(dest_name); % indx in Layers array
                % check if there source has multiple inputs (concat layer, unpooling ...)
                if length(dest) > 1
                    % unpooling option
                    if isa(obj.Layers{dest_indx}, 'MaxUnpooling2DLayer')
                        destP = dest{2};
                        if strcmp(destP, 'in')
                            obj.input_sets{dest_indx} = outSet; % store layer input
                        elseif strcmp(destP, 'indices')
                            obj.Layers{dest_indx}.MaxPoolIndx = obj.Layers{source_indx}.MaxIndx;
                        elseif strcmp(destP, 'size')
                            obj.Layers{dest_indx}.MaxPoolSize = obj.Layers{source_indx}.InputSize;
                        else
                            error("Destination not valid ("+string(obj.Connections.Destination(i))+")");
                        end
                    % concatenation layers (2 or more inputs with specific order {in1, in2, ... inX})
                    elseif contains(class(obj.Layers{dest_indx}), 'Concatenation')
                        destP = dest{2};
                        if startsWith(destP, 'in')
                            input_num = str2double(destP(3:end));
                        else
                            error("Destination not valid ("+string(obj.Connections.Destination(i))+")");
                        end
                        obj.input_sets{dest_indx}{input_num} = outSet;
                    % concatenation layers (2 or more inputs with specific order {in1, in2, ... inX})
                    elseif isa(obj.Layers{dest_indx}, 'DepthConcatenationLayer')
                        destP = dest{2};
                        if startsWith(destP, 'in')
                            input_num = str2double(destP(3:end));
                        else
                            error("Destination not valid ("+string(obj.Connections.Destination(i))+")");
                        end
                        obj.input_sets{dest_indx}{input_num} = outSet;
                    elseif isa(obj.Layers{dest_indx}, 'AdditionLayer')
                        destP = dest{2};
                        if startsWith(destP, 'in')
                            input_num = str2double(destP(3:end));
                        else
                            error("Destination not valid ("+string(obj.Connections.Destination(i))+")");
                        end
                        obj.input_sets{dest_indx}{input_num} = outSet;
                    else
                        error("Destination not valid ("+string(obj.Connections.Destination(i))+")");
                    end
                else
                    obj.input_sets{dest_indx} = outSet;
                end
            end
            % Check if last layer is executed (default is last layer is not executed, but in the 
            % case of the PixelClassificationLayer this is necessary)
            % Assume last layer in array is the output layer
            if isempty(obj.reachSet{end})
                inSet = obj.reachSet{end-1};
                outSet = obj.Layers{end}.reach(inSet, obj.reachMethod, obj.reachOption, obj.relaxFactor, obj.dis_opt, obj.lp_solver);
                obj.reachSet{end} = outSet;
            end
        end
        
        % Ensure precision for layer parameters and inputs is consistent
        % function validate_precision(obj, inSet)
        % 
        % 
        % end
    end % end helper functions
    
    
    methods % semantic segmentation helper functions
    
         % get possible classes of all pixels
        function pc_rs = getPixelClassReachSet(obj, R, rs_id)
            % pc_rs = getPixelClassReachSet(obj, rs_id)
            % @rs_id: reach set id
            % @pc_rs: pixel class reach set
            
            if rs_id < 1 || rs_id > length(R)
                error('Invalid reach set index');
            end
            
            rs = R(rs_id);
            h = rs.height;
            w = rs.width;
            pc_rs = cell(h, w);
            
            for i=1:h
                for j=1:w
                    pc_rs{i, j} = obj.getPixelClasses(R, rs_id, i, j);
                end
            end      
        end
        
        % get all possible classes of a pixels
        function pc = getPixelClasses(~, R, rs_id, h_id, w_id)
            % pc = getPixelClasses(obj, rs_id, h_id, w_id)
            % @rs_id: reach set id
            % @h_id: height index of the pixel
            % @w_id; width index of the pixel
            % @pc: an array of all possible classes of the pixel x(h_id, w_id)
            
            rs = R(rs_id); 
            if h_id < 1 || h_id > rs.height
                error('Invalid height index');
            end
            
            if w_id < 1 || w_id > rs.width
                error('Invalid weight index');
            end
            
            nc = rs.numChannel;
            xmin = zeros(nc,1);
            xmax = zeros(nc,1);
            
            for i=1:nc
                [xmin(i), xmax(i)] = rs.estimateRange(h_id, w_id, i);
            end
            
            max_xmin = max(xmin);
            c = (xmax >= max_xmin);
            pc = find(c);          
        end
        
        % Plot verification result of segmentation task
        function [f1, f2] = plot_segmentation_output_set(obj, verifiedOutput, segImg, figSize)
            %
            % ---- Syntax ----
            % [f1, f2, f3] = plot_segmentation_output_set(obj, verifiedOutput, segImg)
            %
            % ---- Inputs ----
            % @verifiedOutput: verified output of semantic segmentation NN (output of verify_segmentation)
            % @segImgs: ground truth semantic seg images
            % figSize (optional), e.g. [400, 400]
            % 
            % ---- Outputs ----
            % @f1: figure handle of ground truth plot
            % @f2: figure handle of verified output set plot
            
            % Process inputs
            pl_rs = verifiedOutput; % verified output labels
            pl_gr = segImg; % ground truth segmentation image      
            gr_RGB = label2rgb(pl_gr);
            rs_RGB = label2rgb(pl_rs);

            % Get grounf truth info
            gr_unique = unique(pl_gr);
            m1 = length(gr_unique);
            [IND1,in_map1] = rgb2ind(gr_RGB, m1);
            map1 = obj.getColorMap(pl_gr, IND1, in_map1);
            classes1 = obj.Layers{obj.numLayers}.getClasses(gr_unique); 
            
            % Get verified class info
            rs_unique = unique(pl_rs);
            m2 = length(rs_unique);
            [IND2,in_map2] = rgb2ind(rs_RGB, m2);
            map2 = obj.getColorMap(pl_rs, IND2, in_map2);
            classes2 = obj.Layers{obj.numLayers}.getClasses(rs_unique); 
        
            % Set figure size before saving
            if ~exist("figSize", "var")
                figSize = [400 400];
            end
            
            % Create ground truth plot
            f1 = figure;
            ax1 = gca;
            imshow(gr_RGB);
            colormap(ax1,map1);
            cbh1 = colorbar(ax1);
            xtick = 1/(2*m1):1/m1:1;
            cbh1.Ticks = xtick;               
            cbh1.TickLabels = classes1;
            truesize(figSize);
            
            % Create verified output plot
            f2 = figure;
            ax2 = gca;
            imshow(rs_RGB);
            colormap(ax2,map2);
            cbh2 = colorbar(ax2);
            xtick = 1/(2*m2):1/m2:1;
            cbh2.Ticks = xtick;               
            cbh2.TickLabels = classes2;
            truesize(figSize);
            
        end

        % get color map corresponding to the RGB image output of label2rgb()
        function map = getColorMap(~, C, IND, in_map)
            % @C: is the label index matrix
            % @IND:
            % @in_map:
            % @map: return color_map corresponding to the label index in C
            
            % author: Dung Tran
            % date: 4/23/2020
            
            C_unique = unique(C);
            n = length(C_unique);
            [nC, mC] = size(C);
            map = zeros(n, 3);
            for i=1:n
                cat = C_unique(i);
                flag = 0;
                for j=1:nC
                    for k=1:mC
                        if C(j,k) == cat
                            id = IND(j, k) + 1;
                            map(i, :) = in_map(id, :);
                            flag = 1;
                            break;
                        end
                    end
                    if flag == 1
                        break;
                    end
                end
            end
 
        end

    end

    methods (Static)  % semantic segmentation helper functions

        % get paired max pooling layer name
        function maxpooling_layer_name = getPairedMaxPoolingName(Connections, unpooling_layer_name)
            % @unpooling_layer_name: the name of the unmaxpooling layer
            % @maxpooling_layer_name: the name of the paired max pooling layer           
            
            if isempty(Connections)
                error('No connection table');
            end
            
            if ~ischar(unpooling_layer_name)
                error('Invalid unpooling_layer_name');
            else
                dest_name = sprintf("%s/indices", unpooling_layer_name);
            end
            
            n = size(Connections, 1);
            source_name = [];
            for i=1:n                
                if strcmp(Connections.Destination(i), dest_name)
                    source_name = Connections.Source(i);
                    break;
                end
            end
            
            if isempty(source_name)
                error('Unknown destination name');
            end
            
            maxpooling_layer_name = erase(source_name{1}, "/indices");            
        end
    
    end
    
end

