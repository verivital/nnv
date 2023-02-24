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
        reachTime = []; % computation time for each layers
        features = {}; % outputs of each layer in an evaluation
        input_vals = {}; % input values to each layer
        input_sets = {}; % input set values for each layer
        dis_opt = []; % display option = 'display' or []
        lp_solver = 'linprog'; % choose linprog as default LP solver for constructing reachable set user can choose 'glpk' or 'linprog' as an LP solver
        
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
        
        % Evaluation of a NN 
        function y = evaluate(obj, x)
            % Evaluate NN given an input sample
            % y = NN.evaluate(x)
            % @x: input vector x
            % @y: output vector y
            
            % Evaluate layer-by-layer based on connection graph
            for i=1:height(obj.Connections)
                if i > 1
                    x = obj.input_vals{obj.Connections.Source(i)}; % 
                end
                y = obj.Layers{obj.Connections.Source(i)}.evaluate(x);
                if i < height(obj.Connections) % (Last number in destinations is the output)
                    obj.input_vals{obj.Connections.Destination(i)} = y; % store layer input 
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
                if strcmp(obj.dis_opt, 'display')
                    fprintf('\nPerforming analysis for Layer (%s) \n', obj.Layers{obj.Connections.Source(i)});
                end
                rs = obj.reachSet{obj.Connections.Source(i)}; 
                % Compute reachable set layer by layer
                start_time = tic;
                rs_new = obj.Layers{obj.Connections.Source(i)}.reach(rs, obj.reachMethod, obj.reachOption, obj.relaxFactor, obj.dis_opt, obj.lp_solver);
                obj.reachTime(i) = toc(start_time);
                % Store computed reach set
                if i < height(obj.Connections) % (Last number in destinations is the output)
                    obj.reachSet{obj.Connections.Destination(i)} = rs_new; % store layer input 
                end
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
            result = verify_specification(Y, property); 

            % Modify in case of unknown and exact
            if result == 2
                if contains(net.reachMethod, "exact")
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
            [R, ~] = obj.reach(I, reachOptions);        
                      
            % censure correct set format for checking safety
            n = length(R);
            m = length(U);
            R1 = Star;
            if ~isa(R, "Star")
                for i=1:n
                    R1(i) = R(i).toStar; % transform to star sets
                end
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
                        if ~isempty(S) && strcmp(method, 'exact-star')
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
                    end
                end
            end
        end

    end % end verification methods


    methods % specific NN type  methods
        
        % verify robustness of semantic segmentation tasks
        function [riou, rv, rs, n_rb, n_mis, n_unk, n_att, ver_rs] = verify_segmentation(obj, in_images, ground_truths, reachOptions)
            % 
            % --- Syntax ----
            % [riou, rv, rs, n_rb, n_mis, n_unk, n_att, ver_rs] = verify(obj, in_images, ground_truths, reachOptions)
            %
            % ---- Inputs ----
            % @in_images: an array of input set(s)
            % @ground_truths: an array of ground truth images (without attack)
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
            
            % Process inputs
            n1 = length(ground_truths);
            n2 = length(in_images);
            if n1 ~= n2
                error("Inconsistent number of ground truth images and input sets");
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
            
            % compute ground truth output segmentation image
            gr_seg_ims = obj.evaluate_parallel(ground_truths, nCores);
            
            % compute number of correctly classified pixels
            n_pixels = obj.OutputSize(1) * obj.OutputSize(2);
                        
            % obtain verify rearch set
            % 'unknown': the label of the pixel is unknown, may be correct may be not
            % 'misclass': the label of the pixel is misclassified 
            
            % we introduce 2 more classes for the reachable set
            % 1) unknown class
            % 2) misclassification class
                      
            if obj.numCores > 1
                parfor i=1:n1

                    gr_seg_im = gr_seg_ims{i};
                    pc_rs = obj.getPixelClassReachSet(i);
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
                                    ver_im(j,k) = obj.numClasses; % misclass (unrobust pixel)
                                    n_mis_ct = n_mis_ct + 1;
                                end
                            else
                                c = (pc == gr_seg_im(j, k));
                                if sum(c) ~= 0
                                    n_unk_ct = n_unk_ct + 1;
                                    ver_im(j,k) = obj.numClasses - 1; % unknown pixel
                                else
                                    ver_im(j,k) = obj.numClasses; % unrobust pixel
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

                    gr_seg_im = gr_seg_ims{i};
                    pc_rs = obj.getPixelClassReachSet(i);
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
                                    ver_im(j,k) = obj.numClasses; % misclass (unrobust pixel)
                                    n_mis_ct = n_mis_ct + 1;
                                end
                            else
                                c = (pc == gr_seg_im(j, k));
                                if sum(c) ~= 0
                                    n_unk_ct = n_unk_ct + 1;
                                    ver_im(j,k) = obj.numClasses - 1; % unknown pixel
                                else
                                    ver_im(j,k) = obj.numClasses; % unrobust pixel
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
            
            
            obj.verifiedOutputSet = ver_rs;
            obj.groundTruthSegIms = gr_seg_ims;
            obj.RIoU = riou;
            obj.RV = rv;
            obj.RS = rs;
            
            obj.numRbPixels = n_rb;
            obj.numMisPixels = n_mis;
            obj.numUnkPixels = n_unk;
            obj.numPixels = n_pixels;
            obj.numAttPixels = n_att;
            
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
                result = onj.checkRobust(Y(end), target);
            end

        end 

    end % end methods
    

    % helper functions
    methods (Access = protected) % not to be accessible by user
        
        % Check reachability options defined are allowed
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
            g = zeros(height(G),1);

            % Create HalfSapce to define robustness specification
            Hs = HalfSpace(G, g);
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
        
    end % end helper functions
    
    % semantic segmentation helper functions
    methods (Access = protected) 
    
         % get possible classes of all pixels
        function pc_rs = getPixelClassReachSet(obj, rs_id)
            % pc_rs = getPixelClassReachSet(obj, rs_id)
            % @rs_id: reach set id
            % @pc_rs: pixel class reach set
            
            if isempty(obj.reachSet)
                error('Empty reach set, please do reachability analysis first');
            end
            
            if rs_id < 1 || rs_id > length(obj.reachSet)
                error('Invalid reach set index');
            end
            
            rs = obj.reachSet(rs_id);
            h = rs.height;
            w = rs.width;
            pc_rs = cell(h, w);
            
            for i=1:h
                for j=1:w
                    pc_rs{i, j} = obj.getPixelClasses(rs_id, i,j);
                end
            end      
        end
        
        % get all possible classes of a pixels
        function pc = getPixelClasses(obj, rs_id, h_id, w_id)
            % pc = getPixelClasses(obj, rs_id, h_id, w_id)
            % @rs_id: reach set id
            % @h_id: height index of the pixel
            % @w_id; width index of the pixel
            % @pc: an array of all possible classes of the pixel x(h_id, w_id)

            if isempty(obj.reachSet)
                error('Reach Set is empty, please perform reachability first');
            end
            
            rs = obj.reachSet(rs_id); 
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

    end
    
end

