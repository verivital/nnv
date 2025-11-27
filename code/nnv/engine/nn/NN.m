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
        % lp_solver = 'linprog'; % choose linprog as default LP solver for constructing reachable set user can choose 'glpk' or 'linprog' as an LP solver
        lp_solver = 'gurobi';
        matlabnet = []; % the matlab network this NN was created from
        sampling_workers = 5;
        large_set_threshold_sec = 50;   % if a set takes longer than this
            % time to be processed by any layer, declare it to be a "large
            % set", consequently activating adjustment of parallel workers
            % for subsequent layers.
        free_mem_B_before_verify_specification = [];
        free_swap_B_before_verify_specification = [];
        
        % To facilitate graph computation flow
        name2indx = []; % Match name to index in nnvLayers list
        poly_get_mem_using_numel_C = [];  % polynomial to get memory usage
            % of a linear program using number of elements of constraint
            % matrix C. Evaluate polynomial using polyvaln.
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
            obj.fit_poly_to_mem_consumption_of_LP;
                      
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

            % Change if want to execute on GPU
            if isfield(reachOptions, 'device')
                if strcmp(reachOptions.device, 'gpu')
                    obj = obj.params2gpu;
                    inputSet = inputSet.changeDevice('gpu');
                end
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
                else
                    obj.dis_opt = [];
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

            % ensure NN parameters and input set share same precision
            inputSet = obj.consistentPrecision(inputSet); % change only input, this can be changed in the future
            
            % Perform reachability based on connections or assume no skip/sparse connections
            if isempty(obj.Connections)
                outputSet = obj.reach_noConns(inputSet);
            else 
                outputSet = obj.reach_withConns(inputSet, '', reachOptions = reachOptions);
            end

        end
        
        function outputSet = reachProb_ImageStar(obj, IS, reachOptions)


            pe = pyenv;
            py_dir = pe.Executable;


            if isa(IS, 'ImageStar')
                if isempty(IS.im_lb)
                    [LB, UB] = getRanges(IS);
                else
                    LB = IS.im_lb;
                    UB = IS.im_ub;
                end

            elseif isa(IS, 'Star')
                if isempty(IS.state_lb)
                    [LB, UB] = getBox(IS);
                else
                    LB = IS.state_lb;
                    UB = IS.state_ub;
                end

            else
                error('The input must be a Star or Image_Star object.');
            end

            if isfield(reachOptions, 'coverage')
                coverage = reachOptions.coverage;
            else
                coverage = 0.99;
            end
            if isfield(reachOptions, 'confidence')
                confidence = reachOptions.confidence;
            else
                confidence = 0.99;
            end

            if isfield(reachOptions, 'device')
                train_device = reachOptions.device;
            else
                train_device = 'gpu';
            end
            if isfield(reachOptions, 'epochs')
                train_epochs = reachOptions.epochs;
            else
                train_epochs = 50;
            end

            if isfield(reachOptions, 'train_lr')
                train_lr = reachOptions.train_lr;
            else
                train_lr = 0.01;
            end

            [N_dir , N , Ns] = CP_specification(coverage, confidence, numel(IS.im_lb) , train_device, 'single');


            SizeIn = size(LB);
            SizeOut = size(evaluate(obj, LB));
            height = SizeIn(1);
            width = SizeIn(2);


            if isfield(reachOptions, 'indices')
                indices = reachOptions.indices;
            else
                [J,I] = ndgrid(1:width,1:height);
                indices = [I(:), J(:)];
            end


            if isfield(reachOptions, 'mode')
                train_mode = reachOptions.mode;
            else
                train_mode = 'Linear';
            end

            if isfield(reachOptions, 'surrogate_dim')
                surrogate_dim = reachOptions.surrogate_dim;
            else
                surrogate_dim = [-1, -1];
            end

            if isfield(reachOptions, 'threshold_normal')
                threshold_normal = reachOptions.threshold_normal;
            else
                threshold_normal = 1e-5;
            end

            params = struct;
            params.epochs = train_epochs;
            params.lr = train_lr;
            params.trn_batch = floor(N_dir/3);
            params.dims = surrogate_dim;
            params.N_dir = N_dir;
            params.Nt = N;
            params.Ns = Ns;
            % params.files_dir = files_dir;
            params.threshold_normal = threshold_normal;
            params.guarantee = coverage;
            params.py_dir = py_dir;
            params.inputFormat = 'default';

            The_class = ProbReach_ImageStar(obj, LB, UB, indices, SizeOut, train_mode, params);

            outputSet = The_class.ProbReach();

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
            
            outputSet = obj.reach_withConns(inputSet, 'sequence');
            
        end
    
    end
    
    
    methods % secondary methods (verification, safety, robustness...)
        
        % Verify a VNN-LIB specification
        function [result, X] = verify_vnnlib(obj, propertyFile, reachOptions, needReshape)
            
            arguments
                obj 
                propertyFile 
                reachOptions 
                needReshape = 0;
            end
            
            % Load specification to verify
            property = load_vnnlib(propertyFile);
            lb = property.lb;
            ub = property.ub;
            
            if isfield(reachOptions, 'single_average_input') && reachOptions.single_average_input
                if ~isa(lb, 'cell')
                    lb = (lb + ub)/2;
                    ub = lb;
                else
                    lb = (lb{1} + ub{1})/2;
                    ub = lb;
                end
            end
            
            first_layer = obj.Layers(1, 1);
            inputSize = first_layer{1}.InputSize;
            
            % Format bounds into correct dimensions
            % (using code from run_vnncomp2024_instance.m)
            if needReshape == 1
                lb = permute(lb, [2 1 3]);
                ub = permute(ub, [2 1 3]);
            elseif needReshape == 2
                newSize = [inputSize(2) inputSize(1) inputSize(3:end)];
                lb = reshape(lb, newSize);
                lb = permute(lb, [2 1 3 4]);
                ub = reshape(ub, newSize);
                ub = permute(ub, [2 1 3 4]);
            end
            
            % Create reachability parameters and options
            if contains(reachOptions.reachMethod, "zono")
                X = ImageZono(lb, ub);
            else
                X = ImageStar(lb, ub);
            end
            
            % Compute reachability
            if isfield(reachOptions, 'no_reachability_use_given_output_set') && ~isempty(reachOptions.no_reachability_use_given_output_set)
                Y = reachOptions.no_reachability_use_given_output_set;
            else
                Y = obj.reach(X, reachOptions); % Seems to be working
            end
            result = verify_specification(Y, property.prop, reachOptions); 

            % Modify in case of unknown and exact
            if result == 2
                if contains(obj.reachMethod, "exact")
                    result = 0;
                end
            end
    
        end
        
        % Check robustness of output set given a target label
        function rb = checkRobust(obj, outputSet, target, reachOptions)
            % rb = checkRobust(~, outputSet, target)
            % 
            % @outputSet: the outputSet we need to check
            %
            % @target: the correct_id of the classified output or a halfspace defining the robustness specification
            % @rb: = 1 -> robust
            %      = 0 -> not robust
            %      = 2 -> unknown
            
            arguments
                obj 
                outputSet 
                target 
                reachOptions = []
            end
            
            % Process set
            if ~isa(outputSet, "Star")
                nr = length(outputSet);
                R = Star;
                for s=1:nr
                    R(s) = outputSet(s).toStar;
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
            reachOptions.free_mem_B_before_verify_specification = obj.free_mem_B_before_verify_specification;
            reachOptions.free_swap_B_before_verify_specification = obj.free_swap_B_before_verify_specification;
            rb = verify_specification(R, target, reachOptions);
            % if rb ~= verify_specification_old(R, target, reachOptions)
            %     error("Mismatch of robustness verification between old and new code!")
            % end

        end % end check robust
        
        % Verify robutness of a NN given an input set and target (id or HalfSpace)
        function result = verify_robustness(obj, inputSet, reachOptions, target)
            % Compute reachable set
            R = obj.reach(inputSet, reachOptions);
            % Check robustness
            result = obj.checkRobust(R, target, reachOptions);
            % Modify in case of unknown and exact
            if result == 2
                if contains(obj.reachMethod, "exact")
                    result = 0;
                end
            end
            if isfield(reachOptions, 'disp_intersection_result')
                if isfield(reachOptions, 'dis_opt') && strcmp(reachOptions.dis_opt, 'display')
                    fprintf('\n\n')
                end
                fprintf("%s \t ", datetime('now'));
                if result == 1
                    disp('Neural network is ROBUST by intersection.')
                else
                    disp('Neural network robustness UNKNOWN by intersection.')
                end
            end
        end
        
        % Verify robustness of a neural network given a single input
        function [robustness_result_by_intersection, lb, ub, compare_bounds, mismatching_bounds] = verify_robustness_for_3dim_img(obj, reachOptions, args)
            
            arguments
                obj
                reachOptions
                args.input
                args.target = []
                args.matlabnet = []
                args.correct_output = []
                args.plot = false
                args.compare_with_sampling = false
                args.tolerance_for_comparison_with_sampling_bounds = 1e-4
                args.samples_per_pert = 10
                args.netbs = []
                % args.numWorkers_for_par_sampl = net.sampling_workers
            end
            
            if isempty(args.matlabnet)
                args.matlabnet = obj.matlabnet;
            end
            target = args.target;
            [~, pred] = max(predict(args.matlabnet, args.input));
            if ~isempty(target) && pred ~= target
                error('Predicted class different from supplied target class!')
            end
            target = pred;
            
            V = args.input;
            V(:,:,:,2) = zeros(size(args.input));
            
            C = 0;
            d = 0;
            pred_lb = -1;
            pred_ub = 1;
            IS = ImageStar(V, C, d, pred_lb, pred_ub);
            
            robustness_result_by_intersection = obj.verify_robustness(IS, reachOptions, target);
            
            if nargout > 1 % && args.compare_with_sampling
                [lb, ub, robustness_result_by_overlap] = obj.get_ranges_and_plot( ...
                    correct_output = args.correct_output, ...
                    plot = args.plot, ...
                    target = target);
                
                if isfield(reachOptions, 'disp_overlap_result')
                    if robustness_result_by_overlap == 1
                        disp('Neural network is ROBUST by overlap.')
                    else
                        disp('Neural network robustness UNKNOWN by overlap.')
                    end
                end
                
                if args.compare_with_sampling
                    [compare_bounds, mismatching_bounds] = obj.sample_weight_perturbed_nns(input = args.input, ...
                                                    samples_per_pert = args.samples_per_pert, ...
                                                    lb_out = lb, ...
                                                    ub_out = ub, ...
                                                    robustness_result = robustness_result_by_intersection, ...
                                                    matlabnet = args.matlabnet, ...
                                                    netbs = args.netbs, ...
                                                    tolerance_for_comparison_of_bounds = args.tolerance_for_comparison_with_sampling_bounds);
                                                    % , ...
                                                    % numWorkers_for_par_sampl = args.numWorkers_for_par_sampl);
                end
            end
        end
        
        function robustness_result = is_robust_by_overlap(obj, lb, ub, target)
            ub(target) = -Inf;
            if any(lb(target) < ub)
                robustness_result = false;
            else
                robustness_result = true;
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
    

    methods     % helper functions
        
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
        
        % Ensure input and parameter precision is the same
        function inputSet = consistentPrecision(obj, inputSet)
            % (assume parameters have same precision across layers)
            % approach: change input precision based on network parameters
            inputPrecision = class(inputSet(1).V);
            netPrecision = 'double'; % default
            for i=1:length(obj.Layers)
                if isa(obj.Layers{i}, "FullyConnectedLayer") || isa(obj.Layers{i}, "Conv2DLayer")
                    netPrecision = class(obj.Layers{i}.Weights);
                    break;
                end
            end
            if ~strcmp(inputPrecision, netPrecision)
                % input and parameter precision does not match
                warning("Changing input set precision to "+string(netPrecision));
                for i = 1:length(inputSet)
                    inputSet(i) = inputSet(i).changeVarsPrecision(netPrecision);
                end
            end
        end

        % Change paramters to gpu
        function obj = params2gpu(obj)
            % change the parameters layer by layer
            for i = 1:length(obj.Layers)
                gpuLayer = obj.Layers{i}.toGPU;
                obj.Layers{i} = gpuLayer;
            end
        end

        % Change params precision
        function obj = changeParamsPrecision(obj, precision)
            % change the parameters layer by layer
            for i = 1:length(obj.Layers)
                pLayer = obj.Layers{i}.changeParamsPrecision(precision);
                obj.Layers{i} = pLayer;
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
        function start_pool(obj, args)
            
            arguments
                obj 
                args.allow_for_1_worker = false
            end

            if obj.numCores > 1 || ((obj.numCores == 1) && args.allow_for_1_worker)
                poolobj = gcp('nocreate'); % If no pool, do not create new one.
                if isempty(poolobj)
                    parpool('local', obj.numCores); 
                else
                    if poolobj.NumWorkers ~= obj.numCores
                        delete(poolobj); % delete the old poolobj
                        if obj.numCores > 1 || ((obj.numCores == 1) && args.allow_for_1_worker)
                            parpool('local', obj.numCores); % start the new one with new number of cores
                        end
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
        function outSet = reach_withConns(obj, inSet, reachType, args)
            
            arguments
                obj 
                inSet 
                reachType = 'default'
                args.reachOptions = struct
            end
            
            reachOptions = args.reachOptions;
            
            % Initialize variables to store reachable sets and computation time
            obj.reachSet = cell(1, obj.numLayers);
            obj.reachTime = zeros(1, obj.numLayers);
            obj.input_sets = cell(1, height(obj.Connections)); % store input reach sets for each layer
            if strcmp(obj.dis_opt, 'display')
                fprintf('\nPerform reachability analysis for the network %s...', obj.Name);
            end
            
            high_gpu_mem_in_last_layer = 0;
            if isfield(reachOptions, 'layer_specific_numCores')
                obj.adjust_numCores("start", reachOptions);
            end
            too_big_for_gpu = 0;
            
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
                    if isfield(reachOptions, 'delete_old_sets') && reachOptions.delete_old_sets	  % delete reachable sets of preceding layers to save memory when processing huge star sets
                       obj.input_sets{source_indx} = {};   % deleting set may cause issues
                    end
                end
                
                % 3) reach layer
                exec_len = length(obj.reachSet);
                % ensure layer has not been eexcuted yet
                if exec_len >= source_indx && isempty(obj.reachSet{source_indx})
                    if strcmp(reachOptions.reachMethod, 'approx-star')
                        inSet1 = inSet(1);
                        if isa(inSet1, 'cell')
                            inSet1 = inSet1{1};
                        end
                        if isfield(reachOptions, 'layer_specific_numCores')
                            obj.adjust_numCores(class(obj.Layers{source_indx}), reachOptions, size(inSet1.C), size(inSet1.V), underlyingType(inSet1.V));
                        end
                    end
                    
                    if strcmp(obj.dis_opt, 'display')
                        fprintf("\n\n%s \t ", datetime('now'));
                        fprintf('Performing analysis for Layer %d (%s) with %d input set(s)... ', i-1, source, length(inSet));
                        if isscalar(inSet)
                            fprintf("\n\t Dimensions: %d. ", inSet(1).numPred);
                            fprintf("Number of constraints: %d. ", size(inSet(1).C, 1));
                            fprintf("Product: %d. ", numel(inSet(1).C));
                            poolobj = gcp('nocreate'); % If no pool, do not create new one.
                            if isempty(poolobj)
                                numParWorkers = 1;
                            else
                                numParWorkers = poolobj.NumWorkers;
                            end
                            fprintf("Parallel workers: %d. ", numParWorkers);
                            disp('');
                        end
                    end
                    
                    switch_devices_for_efficiency = (isa(obj.Layers{source_indx}, 'MaxPooling2DLayer') ...
                                                  || isa(obj.Layers{source_indx}, 'ReluLayer') ...
                                                  || isa(obj.Layers{source_indx}, 'DepthConcatenationLayer') ...
                                                  ) && isfield(reachOptions, 'device') ...
                                                    && strcmp(reachOptions.device, 'gpu');
                    
                    if switch_devices_for_efficiency
                        ImageStar.change_device_multiple_sets(inSet, 'cpu');
                    end
                    
                    if too_big_for_gpu && isa(inSet, 'cell')
                        for k = 1:length(inSet)
                            sets = inSet{k};
                            for m = 1:length(sets)
                                if isa(sets(m).V, 'gpuArray')
                                    sets(m).changeDevice('cpu');
                                end
                            end
                        end
                    end
                    
                    if isfield(reachOptions, "disp_layer_type")
                        fprintf("\n\t\t Layer type: %s. Input: ", class(obj.Layers{source_indx}));
                        if ~isa(inSet, 'ImageStar')
                            for k = 1:length(inSet)
                                disp(inSet{k})
                                if isempty(inSet{k}.V)
                                    error("Input set " + k + " is empty!")
                                end
                            end
                        else
                            disp(inSet)
                            if isempty(inSet.V)
                                error("Input set is empty!")
                            end
                        end
                    end
                    
                    t = tic;
                    if strcmp(reachType, 'sequence')
                        outSet = obj.Layers{source_indx}.reachSequence(inSet, obj.reachMethod, obj.reachOption, obj.relaxFactor, obj.dis_opt, obj.lp_solver);
                    else
                        outSet = obj.Layers{source_indx}.reach(inSet, obj.reachMethod, obj.reachOption, obj.relaxFactor, obj.dis_opt, obj.lp_solver);
                    end
                    
                    obj.reachTime(source_indx) = toc(t);
                    if strcmp(obj.dis_opt, 'display')
                        fprintf(' Done in %g s.', obj.reachTime(source_indx));
                    end
                    
                    if isa(outSet(1).V, 'gpuArray')
                        if NN.get_nvidia_gpu_used_memory_frac > 0.5
                            poolobj = gcp('nocreate');
                            if ~isempty(poolobj)
                                n = poolobj.NumWorkers;
                                delete(gcp('nocreate'));
                                parpool('local', n);
                            end
                            if ~high_gpu_mem_in_last_layer
                                high_gpu_mem_in_last_layer = 1;
                            else
                                outSet(1).too_big_for_gpu = 1;
                                too_big_for_gpu = 1;
                                ImageStar.change_device_multiple_sets(outSet, 'cpu');
                            end
                        else
                            high_gpu_mem_in_last_layer = 0;
                        end
                    end
                    
                    for k = 1:length(outSet)
                        if outSet(k).too_big_for_gpu
                            obj = obj.params2cpu;
                            too_big_for_gpu = 1;
                            switch_devices_for_efficiency = 0;
                            break;
                        end
                    end
                    
                    if switch_devices_for_efficiency
                        ImageStar.change_device_multiple_sets(outSet, 'gpu');
                    end
                    
                    if (~isfield(reachOptions, 'delete_old_sets')) || (~reachOptions.delete_old_sets) || (source_indx > length(obj.Layers) - 2)
                       obj.reachSet{source_indx} = outSet;
                    else
                       obj.reachSet{source_indx} = ["Fake reachable set, to save memory."];
                    end
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
                if strcmp(reachType, 'sequence')
                    outSet = obj.Layers{end}.reachSequence(inSet, obj.reachMethod, obj.reachOption, obj.relaxFactor, obj.dis_opt, obj.lp_solver);
                else
                    outSet = obj.Layers{end}.reach(inSet, obj.reachMethod, obj.reachOption, obj.relaxFactor, obj.dis_opt, obj.lp_solver);
                end
                obj.reachSet{end} = outSet;
            end
            
            if isfield(reachOptions, 'layer_specific_numCores')
                [free_mem_B, mem_margin_B, mem_per_LP_B] = obj.adjust_numCores("end", reachOptions, size(outSet(1).C) + [1 0]);
                if isfield(reachOptions, 'dis_opt') && strcmp(reachOptions.dis_opt, 'display')
                    fprintf("\n\n%s \t ", datetime('now'));
                    fprintf("Free memory: %.4g GB, ", free_mem_B/2^30)
                    fprintf("Margin: %.4g GB, ", mem_margin_B/2^30)
                    fprintf("Memory per LP: %.4g GB\n", mem_per_LP_B/2^30)
                    obj.free_mem_B_before_verify_specification = free_mem_B;
                    obj.free_swap_B_before_verify_specification = NN.get_free_swap_B;
                end
            end
        end
        
        % return the indices of layers of the specified types
        function layers = get_layer_indices(obj, layer_types)
            layers = [];
            for l = 1:length(obj.Layers)
                for n = 1:length(layer_types)
                    if isa(obj.Layers{l}, layer_types(n))
                        layers = [layers l];
                    end
                end
            end
        end
        
        % sample a neural network's perturbation space by plugging the
        % perturbed weights into the network and testing the accuracy over
        % the supplied input. It is assumed that the supplied input is
        % classified correctly by the unperturbed neural network, and
        % predictions different from that of the unperturbed neural network
        % are considered misclassifications.
        function [compare_bounds, mismatching_bounds] = sample_weight_perturbed_nns(obj, args)
            
            arguments
                obj                             % the NN object
                args.input                      % the input to the NN
                args.lb_out                     % lower bounds on the
                    % output of the neural network obtained using NNV
                args.ub_out                     % upper bounds on the
                    % output of the neural network obtained using NNV
                args.robustness_result          % robustness result
                    % obtained using NNV
                args.samples_per_pert = 10      % samples per perturbation
                args.display = 1                % whether to display results
                args.netbs = []                 % bytestream based copy of 
                    % the NNV NN object
                % args.numWorkers_for_par_sampl = 5
                args.tolerance_for_comparison_of_bounds = 1e-4
            end
            
            samples_per_pert = args.samples_per_pert;    % no. of samples per weight
            
            % gather the weight perturbations from the various layers of
            % the neural network
            weightPerturbSpecs = [];
            for l = 1:length(obj.Layers)
                if isa(obj.Layers{l}, 'FullyConnectedLayer') || isa(obj.Layers{l}, 'Conv2DLayer')
                    weightPerturbSpecsThisLayer = obj.Layers{l}.weightPerturb;
                    sz = size(weightPerturbSpecsThisLayer);
                    if sz(2) > 0
                        weightPerturbSpecs = [weightPerturbSpecs; l*ones(sz(1),1) weightPerturbSpecsThisLayer];
                    end
                elseif isprop(obj.Layers{l}, 'weightPerturb')
                    error("Layer " + obj.Layers{l}.Name + " not supported for sampling yet.")
                end
            end
            
            numPerts = size(weightPerturbSpecs, 1);
            weightUB = weightPerturbSpecs(:, end);
            weightLB = weightPerturbSpecs(:, end - 1);
            deltasWeightPerturb = (weightUB - weightLB)/(samples_per_pert - 1);
                % deltas of weight perturbations
            
            powersNSamplesPerPert = samples_per_pert.^(numPerts - 1 :-1: 0).';   % powers of n
            
            correct_output = obj.evaluate(args.input);
            correct_output = correct_output(:);
            [~, target] = max(correct_output);
            
            lb_netOut = ones(size(correct_output))*Inf;
            ub_netOut = -lb_netOut;
            
            incorrectly_classified_samples = 0;
            
            total_combs = samples_per_pert^numPerts;
            % old_numCores = obj.numCores;
            % obj.numCores = args.numWorkers_for_par_sampl;
            % obj.start_pool;
            parfor m = 1:total_combs      % perturbation combination number
            % for m = 1:total_combs      % perturbation combination number
                multDeltasWeightPerturb = mod(ceil(m./powersNSamplesPerPert) - 1, samples_per_pert);  % multiples of deltas of weight perturbations
                weightPerturb = multDeltasWeightPerturb.*deltasWeightPerturb + weightLB;   % weight perturbations
                
                tmp_net = getArrayFromByteStream(args.netbs);
                
                for k = 1 : numPerts
                    l = weightPerturbSpecs(k, 1);
                    ind = weightPerturbSpecs(k, 2);
                    if isa(obj.Layers{l}, 'FullyConnectedLayer') || isa(obj.Layers{l}, 'Conv2DLayer')
                        n = numel(obj.Layers{l}.Weights);
                        if ind <= n     % perturb weights matrix
                            tmp_net.Layers{l}.Weights(ind) = tmp_net.Layers{l}.Weights(ind) + weightPerturb(k);
                        else            % perturb bias matrix
                            ind = ind - n;
                            tmp_net.Layers{l}.Bias(ind) = tmp_net.Layers{l}.Bias(ind) + weightPerturb(k);
                        end
                    else
                        error(['Layer ' class(obj.Layers{l}) ' not supported in sampling!'])
                    end
                end
                
                netOut = tmp_net.evaluate(args.input);
                netOut = netOut(:);
                [~, predThisComb] = max(netOut);
                if predThisComb ~= target
                    incorrectly_classified_samples = incorrectly_classified_samples + 1;
                end
                
                lb_netOut = min(lb_netOut, netOut);
                ub_netOut = max(ub_netOut, netOut);
            end
            
            % debug_disp
            lb_out = args.lb_out;
            ub_out = args.ub_out;
            compare_bounds = [lb_netOut, ub_netOut, lb_netOut - lb_out, ub_out - ub_netOut, lb_out, ub_out];
                % bounds from sampling, overapproximation of bounds
                % from star sets over bounds from sampling, bounds from
                % star sets
            
            percentage_of_incorrect_classification = incorrectly_classified_samples/total_combs*100;
            fprintf('\n\n------------------------- Comparison with Sampling -------------------------\n\n')
            disp(['Samples incorrectly classified: ' num2str(percentage_of_incorrect_classification) ' %'])
            
            if any(lb_netOut < 0, 'all')
                shift = abs(min(lb_netOut, [], 'all'));
            else
                shift = 0;
            end
            shifted_lb_netOut = lb_netOut + shift;
            shifted_ub_netOut = ub_netOut + shift;
            
            zero_thresh = args.tolerance_for_comparison_of_bounds*shifted_lb_netOut(target);
            
            matching_lbs = abs(compare_bounds(:,3)) < zero_thresh;
            matching_ubs = abs(compare_bounds(:,4)) < zero_thresh;
            matching_bounds = matching_lbs & matching_ubs;
            fprintf('Tolerance: %g \t ', args.tolerance_for_comparison_of_bounds)
            mismatching_bounds = [];
            if all(matching_bounds)
                fprintf('Bounds match.\n')
            else
                fprintf('Bounds MISMATCH! Mismatching bounds:\n\n')
                mismatching_bounds = [find(~matching_bounds), compare_bounds(~matching_bounds, :)];
                disp(mismatching_bounds)
            end
            
            if (incorrectly_classified_samples > 0) && (args.robustness_result == 1)
                % one or more samples resulted in incorrect classification, even though
                % neural network was determined to be robust by nnv
                error('Robustness mismatch!')
            end
            
            % net.numCores = old_numCores;
        end
        
        % get output ranges of the reachable set(s) of the neural network
        function [lb, ub, robustness_result_by_overlap] = get_ranges_and_plot(obj, args)
            
            arguments
                obj
                args.correct_output = []
                args.plot = false
                args.target = []
                % args.x_labels = []
                args.lp_solver = obj.lp_solver
            end
            
            robustness_result_by_overlap = nan;

            % Get output reachable set
            R = obj.reachSet{end};
            for k = 1:length(R)
                R(k).changeDevice('cpu');
            end
            
            [lb, ub] = R(1).getRanges(args.lp_solver, 'parallel');
            lb = squeeze(lb);
            ub = squeeze(ub);
            if ~isempty(args.target)
                robustness_result_by_overlap = obj.is_robust_by_overlap(lb, ub, args.target);
            end
            
            % obj.start_pool;
            progress = zeros(length(R), 1);
            progress(1) = 1;
            parfor k = 2:length(R)
                % Get (overapproximate) ranges for each output index
                [lb_out_k, ub_out_k] = R(k).getRanges(args.lp_solver);
                lb_out_k = squeeze(lb_out_k);
                ub_out_k = squeeze(ub_out_k);
                if ~isempty(args.target)
                    robustness_result_by_overlap = robustness_result_by_overlap & obj.is_robust_by_overlap(lb_out_k, ub_out_k, args.target);
                end
                lb = min(lb, lb_out_k);
                ub = max(ub, ub_out_k);
                progress(k) = 1;
            end
            
            if args.plot
                % Get middle point for each output and range sizes
                mid_range = (lb + ub)/2;
                range_size = ub - mid_range;
                
                % Label for x-axis
                x = [1:length(lb)]';
                
                % Visualize set ranges and evaluation points
                hold off;
                errorbar(x, mid_range, range_size, '.');
                xlim([(x(1) - 0.5) (x(end) + 0.5)]);
                if ~isempty(args.correct_output)
                    hold on;
                    scatter(x, args.correct_output, 'x', 'MarkerEdgeColor', 'r');
                    hold off;
                end
            end
            
        end
        
        function fc_layers = verify_has_only_fc_relu_placeholder_layers(obj, l_start, l_end)
            
            arguments
                obj
                l_start 
                l_end = []
            end
            
            if isempty(l_end)
                fc_layers = obj.get_layer_indices("FullyConnectedLayer");
                l_end = fc_layers(end);
            end
            
            fc_layers = [];
            relu_encountered = 0;
            for l = l_start:l_end
                if isa(obj.Layers{l}, "FullyConnectedLayer")
                    fc_layers = [fc_layers l];
                    relu_encountered = 0;
                elseif isa(obj.Layers{l}, "ReluLayer")
                    if relu_encountered
                        error("Two ReLU layers without FullyConnectedLayer between them not supported.")
                    end
                    relu_encountered = 1;
                elseif isa(obj.Layers{l}, "PlaceholderLayer")
                else
                    error_struct.message = ['Unsupported layer ' class(obj.Layers{l}) '.'];
                    error_struct.identifier = "NNV:non_supported_layers_in_fc_relu_only_method";
                    error(error_struct);
                end
            end
            
            if ~isa(obj.Layers{l_start}, "FullyConnectedLayer")
                error(['First layer should be a fully connected layer.'])
            end
            if ~isa(obj.Layers{l_end}, "FullyConnectedLayer")
                error(['Last layer should be a fully connected layer.'])
            end
            
        end
        
        % for comparison with another approach in ModelStar experiments. Can be removed.
        function [lb, ub] = compute_bounds_weng_2020_fc_layers_only(obj, args)
            % call this function after performing reachability analysis
            % using star sets.
            
            arguments
                obj
                args.first_fc_layer
                args.last_fc_layer = []
                args.pert_vec = []
                args.pert_mag = []
                args.pert_norm = Inf
            end
            
            p = args.pert_norm;
            N = args.first_fc_layer;
            K = args.last_fc_layer;
            pert_vec = args.pert_vec;
            % N is the perturbed layer.
            % pert_vec is the perturbation applied to every row of
            % layer N's weight matrix.
            
            fc_layers = obj.verify_has_only_fc_relu_placeholder_layers(N, K);
            
            W = @(k) obj.Layers{fc_layers(k)}.Weights;
            b = @(k) obj.Layers{fc_layers(k)}.Bias;
            
            l = cell(length(fc_layers), 1);
            u = l;
            alphaL = l;
            betaL = l;
            alphaU = l;
            betaU = l;
            delta = l;
            theta = l;
            
            xN = obj.input_sets{N}.V(:,:,:,1);
            xN = xN(:);
            if isempty(pert_vec)
                pert_vec = args.pert_mag*ones(1, length(xN));
            end
            yN = W(1)*xN + b(1);
            pert = pert_vec*abs(xN);
            l{1} = yN - pert;
            u{1} = yN + pert;
            
            e = norm(pert_vec, p);
            q = 1/(1 - 1/p);
            
            % compute bounds layer-wise i.e., 2nd layer, then 3rd layer, and so on
            for K = 2:length(fc_layers)
                numOutputs = size(W(K), 1);
                
                A = cell(length(fc_layers), 1);
                O = A;
                lambda = A;
                omega = A;
                
                m = K - 1;
                
                alphaU{m} = u{m}./(u{m} - l{m});
                alphaL{m} = alphaU{m};
                 betaU{m} = -l{m};
                
                alphaU{m}(u{m} <= 0) = 0;
                alphaL{m}(u{m} <= 0) = 0;
                 betaU{m}(u{m} <= 0) = 0;
                
                alphaU{m}(l{m} >= 0) = 1;
                alphaL{m}(l{m} >= 0) = 1;
                 betaU{m}(l{m} >= 0) = 0;
                
                u{K} = nan(numOutputs, 1);
                l{K} = u{K};
                
                for j = 1:numOutputs
                    ej = zeros(numOutputs, 1);
                    ej(j) = 1;
                    
                    A{K}(j,:) = ej';
                    O{K}(j,:) = ej';
                    
                    for k = K-1 : -1 : 1
                        A{k} = nan(numOutputs, size(W(k+1), 2));
                        O{k} = A{k};
                        lambda{k} = A{k};
                        omega{k} = A{k};
                        delta{k} = A{k}';
                        theta{k} = A{k}';
                        
                        Wkp1 = W(k+1);
                        for i = 1:size(Wkp1, 2)
                            if A{k+1}(j,:) * Wkp1(:,i) >= 0
                                lambda{k}(j,i) = alphaU{k}(i);
                                 delta{k}(i,j) =  betaU{k}(i);
                            else
                                lambda{k}(j,i) = alphaL{k}(i);
                                 delta{k}(i,j) = 0;     % betaL is always 0
                            end
                            if O{k+1}(j,:) * Wkp1(:,i) >= 0
                                 omega{k}(j,i) = alphaL{k}(i);
                                 theta{k}(i,j) = 0;     % betaL is always 0
                            else
                                 omega{k}(j,i) = alphaU{k}(i);
                                 theta{k}(i,j) =  betaU{k}(i);
                            end
                        end
                        A{k}(j,:) = (A{k+1}(j,:) * W(k+1)).*lambda{k}(j,:);
                        O{k}(j,:) = (O{k+1}(j,:) * W(k+1)).* omega{k}(j,:);
                    end
                    
                    t4 = A{K}(j,:) * b(K);    % delta{K} is always zero
                    t5 = O{K}(j,:) * b(K);    % theta{K} is always zero
                    for k = 1:K-1
                        t4 = t4 + A{k}(j,:) * (b(k) + delta{k}(:,j));
                        t5 = t5 + O{k}(j,:) * (b(k) + theta{k}(:,j));
                    end
                    
                    u{K}(j) =  e * norm(A{1}(j,:), 1) * norm(xN, q) + A{1}(j,:) * W(1)*xN + t4;
                    l{K}(j) = -e * norm(O{1}(j,:), 1) * norm(xN, q) + O{1}(j,:) * W(1)*xN + t5;
                end
            end
            
            lb = l{end};
            ub = u{end};
        end
        
        % for comparison with another approach in ModelStar experiments. Can be removed.
        function image_verified = formal_robust_verify_fc_layers_only(obj, args)
            
            arguments
                obj
                args.first_fc_layer
                args.pert_mag
                args.correct_output
                args.last_fc_layer = []
            end
            
            p = args.pert_mag;
            
            Y_outputs = args.correct_output;
            [~, yPred] = max(Y_outputs);
            
            N = args.first_fc_layer;
            K = args.last_fc_layer;
            fc_layers = obj.verify_has_only_fc_relu_placeholder_layers(N, K);
            K = fc_layers(end);
            
            xN = obj.input_sets{N}.V(:,:,:,1);
            xN = xN(:);     % input to perturbed layer
            
            W = obj.Layers{K}.Weights;      % last layer's weight matrix
            
            image_verified = 1;
            for j = 1:length(Y_outputs)
                if j ~= yPred
                    t1 = norm(W(yPred,:) - W(j,:), 1);
                    if N == K
                        t1 = 2;
                    end
                    
                    t2 = norm(xN, 1);                   % input to the only perturbed layer
                    
                    t3 = 1;
                    for l = fc_layers(2:end - 1)   % all fully-connected layers between the perturbed layer and the last layer
                        t3 = t3*max(norms(obj.Layers{l}.Weights, 1));
                    end
                    
                    ub_on_decrease_in_margin = p*t1*t2*t3;
                    original_margin = Y_outputs(yPred) - Y_outputs(j);
                    if ub_on_decrease_in_margin >= original_margin
                        image_verified = 0;
                        break
                    end
                end
            end
        end
        
        function the_range = get_weights_range(obj, layer_no)
            the_range = range(obj.Layers{layer_no}.Weights, 'all');
        end
        
        function print_layers_info(obj)
            nn_info = analyzeNetwork(obj.matlabnet);
            for l_no = 1:length(obj.Layers)
                l = obj.Layers{l_no};
                fprintf("Layer %2d: %-30s %-25s ", l_no, class(l), l.Name)
                try
                    fprintf("num_Weights: %-10d ", numel(l.Weights));
                    tmp = nn_info.LayerInfo(l_no, 3);
                    tmp = tmp{1, 1};
                    n_neurons = prod(cell2mat(tmp(1)));
                    fprintf("num_Neurons: %-10d", n_neurons);
                catch ME
                    if ~strcmp(ME.identifier, 'MATLAB:noSuchMethodOrField')
                        rethrow(ME)
                    end
                end
                disp(' ')
            end
        end
        
        function [actual_free_mem_B, mem_margin_B, mem_per_LP_B] = adjust_numCores(obj, layer_type, reachOptions, sz_C, sz_V, type_V)
            
            arguments
                obj 
                layer_type 
                reachOptions 
                sz_C = []
                sz_V = 0
                type_V = 'single'
            end
            
            mem_margin = 0.1;   % ignore this percentage of free memory
            actual_free_mem_B = NN.get_free_mem_B;
            mem_V_B = NN.get_required_mem_B(sz_V, type_V);
            mem_margin_B = mem_margin*actual_free_mem_B;
            free_mem_B = actual_free_mem_B - mem_margin_B;
            if ~strcmp(layer_type, "end")
                free_mem_B = free_mem_B - mem_V_B;
            end
            if strcmp(layer_type, 'ReluLayer')
                free_mem_B = free_mem_B - 2*mem_V_B;
            end
            mem_margin_B = actual_free_mem_B - free_mem_B;
            if ~isempty(sz_C)
                mem_per_LP_B = obj.estimate_LP_mem_usage_B(sz_C);
            else
                mem_per_LP_B = 1;
            end
            if free_mem_B < 1
                % error("Insufficient memory to hold the V of the new set with mem_margin " + mem_margin);
                free_mem_B = mem_per_LP_B;    % debug, allow swap usage
            end
            
            if isfield(reachOptions, 'layer_specific_numCores') && reachOptions.layer_specific_numCores.isKey(layer_type)
            
                numParWorkers = reachOptions.layer_specific_numCores(layer_type);
                if numParWorkers <= 0
                    maxWorkers = parcluster('local').NumWorkers;
                    numParWorkers = max(1, min(floor(free_mem_B/mem_per_LP_B), maxWorkers));
                end
                
                old_numCores = obj.numCores;
                obj.numCores = numParWorkers;
                obj.start_pool(allow_for_1_worker = true);
                obj.numCores = old_numCores;
            end
        end
        
        function p = fit_poly_to_mem_consumption_of_LP(obj)
            % data consists of columns of C, rows of C, number of elements
            % of C and memory consumed by one linear program using that C.
            data = [0,      0,      0    ,  0  ;
                    3985,   11851,  0.81 ,  0  ;
                    169256, 256,    3.78 ,  0  ;
                    19752,  3958,   4.98 ,  0  ;
                    19709,  3829,   5.84 ,  0  ;
                    101389, 4025,   6    ,  0  ;
                    147853, 423,    6.38 ,  8  ;
                    19907,  4423,   6.67 ,  17 ;
                    102759, 8135,   15.5 ,  0  ;
                    17797,  53287,  19.6 ,  0  ;
                    16875,  75217,  23   ,  0  ;
                    105288, 15724,  41   ,  4  ;
                    % 105288, 15723,  36.8 ,  0  ;
                    % 214358, 135307, 830  ,  0  ;
                    ];
            original_mem_consump = data(:, 3);
            data(:, 3) = data(:, 3) + data(:, 4);   % bias for fitting to quadratic
            data(:, 3) = data(:, 3)*2^30;
            x = data(:, [1 2]); % input is number of columns and rows of C
            y = data(:, 3); % output is memory consumed per LP
            % p = polyfitn(x,y,1);
            p = polyfitn(x,y,2);
            % [original_mem_consump, polyvaln(p, x)/2^30]   % debug
            
            obj.poly_get_mem_using_numel_C = p;
        end
        
        function mem_B = estimate_LP_mem_usage_B(obj, sz_C)
            % n = [sz_C(2), sz_C(1), prod(sz_C)];
            n = [sz_C(2), sz_C(1)];
            mem_B = max(1, polyvaln(obj.poly_get_mem_using_numel_C, n));
        end
        
        % for ModelStar experiments. Can be removed.
        function net = acas_combine_affine_layer_with_fc(obj, args)
            
            arguments
                obj
                args.delete_affine_layers = 0
            end
            
            wrong_network = 0;
            % wrong_network = wrong_network || ~(isempty(obj.reachSet) && isempty(obj.input_sets));
            
            wrong_network = wrong_network || ...
                            ~isa(obj.Layers{1},  'ImageInputLayer') || ~strcmp(obj.Layers{1}.Normalization, 'none') || ...
                            ~isa(obj.Layers{2}, 'PlaceholderLayer') || ...
                                 obj.Layers{3}.DoOffset || ...
                            ~isa(obj.Layers{4},     'FlattenLayer');
            
            acas_fc_layers = 5:3:23;
            wrong_network = wrong_network || norm(find(arrayfun(@(l) isa(l{1}, 'FullyConnectedLayer'), obj.Layers)).' - acas_fc_layers, 1) ~= 0;
            
            acas_elementwise_affine_layers = 3:3:24;
            wrong_network = wrong_network || norm(find(arrayfun(@(l) isa(l{1}, 'ElementwiseAffineLayer'), obj.Layers)).' - acas_elementwise_affine_layers, 1) ~= 0;
            wrong_network = wrong_network || any(arrayfun(@(l) isa(l{1}, 'ElementwiseAffineLayer') && l{1}.DoScale, obj.Layers));
            
            acas_relu_layers = 7:3:22;
            wrong_network = wrong_network || norm(find(arrayfun(@(l) isa(l{1}, 'ReluLayer'), obj.Layers)).' - acas_relu_layers, 1) ~= 0;
            
            wrong_network = wrong_network || ~isa(obj.Layers{25}, 'PlaceholderLayer') || length(obj.Layers) > 25;
            
            if wrong_network
                error("Network did not satisfy the expected ACAS network structure.")
            end
            
            for fc_layer = acas_fc_layers
                if (fc_layer < acas_fc_layers(end) && any(size(obj.Layers{fc_layer + 1}.Offset) ~= [50 1])) || (fc_layer == acas_fc_layers(end) && any(size(obj.Layers{fc_layer + 1}.Offset) ~= [5 1]))
                    error("Wrong size of offset.")
                end
                obj.Layers{fc_layer}.Bias = obj.Layers{fc_layer}.Bias + obj.Layers{fc_layer + 1}.Offset;
                obj.Layers{fc_layer + 1}.Offset = zeros(size(obj.Layers{fc_layer + 1}.Offset));
                obj.Layers{fc_layer + 1}.DoOffset = 0;
            end
            
            if any(arrayfun(@(l) isa(l{1}, 'ElementwiseAffineLayer') && (l{1}.DoScale || l{1}.DoOffset), obj.Layers))
                error("There is some ElementwiseAffineLayer doing some scaling or offsetting.")
            end
            
            if args.delete_affine_layers
                if ~isempty(obj.reachSet) || ~isempty(obj.input_sets)
                    if any(arrayfun(@(l) obj.input_sets{l} ~= obj.reachSet{l}, acas_elementwise_affine_layers))
                        error("ElementwiseAffineLayer's input and output should have been identical i.e., it should have been combined with the preceding fully connected layer before reachability analysis.")
                    end
                end
                obj.Layers(acas_elementwise_affine_layers) = [];
                if ~isempty(obj.reachSet) || ~isempty(obj.input_sets)
                    obj.reachSet(acas_elementwise_affine_layers) = [];
                    obj.input_sets(acas_elementwise_affine_layers) = [];
                end
            end
            
            net = obj;
            
        end
        
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
    
    methods (Static)    % weight perturbation helper functions
        
        % specify a weight perturbation to be added to a layer's
        % weightPerturb specification matrix
        function add_pert_call_this_function_from_layer(layer, pert_indices, lb, ub)
            % @layer: layer to specify perturbation for
            % @pert_indices: indices of the weight to be perturbed
            % @lb: lower bound of the perturbation
            % @ub: upper bound of the perturbation
            
            sz = size(layer.Weights);
            if length(pert_indices) ~= length(sz)
                error(['Must specify location of perturbation in layer ' layer.Name ' as indices of its weights matrix which has size ' sz])
            end
            
            pert_indices = num2cell(pert_indices);
            linear_ind = sub2ind(sz, pert_indices{:});
            
            layer.weightPerturb = [layer.weightPerturb; linear_ind lb ub];
        end
        
        % specify a perturbation to be applied to the whole layer
        function perturb_whole_layer_call_this_function_from_layer(layer, lb, ub)
            % @lb: lower bound of the perturbation
            % @ub: upper bound of the perturbation
            
            n = numel(layer.Weights);
            
            temp_pert = nan(n, 3);
            temp_pert(:, 1) = [1:n]';
            temp_pert(:, 2) = lb;
            temp_pert(:, 3) = ub;
            
            layer.weightPerturb = [layer.weightPerturb; temp_pert];
        end
        
        % specify a symmetric perturbation to be applied to the whole layer
        % as a fraction of the range of the weights in the layer
        function pert_whole_layer_given_fraction_of_weights_range_call_fromLayer(layer, frac)
            % @frac: fraction of the weights' range to be applied as
            % perturbation to the whole layer
            
            p = frac*range(layer.Weights, 'all');
            NN.perturb_whole_layer_call_this_function_from_layer(layer, -p, p);
        end
        
        function mem = get_free_mem_B()
            [~, mem_str] = system("vmstat -s -S K | grep free\ memory | awk '{print $1}'");
            mem = str2double(mem_str)*2^10;
        end
        
        function idle_cpu = get_idle_cpu()
            [~, idle_cpu] = system("mpstat 1 1 | sed -n '4p' | awk '{print $NF}'");
            idle_cpu = str2double(idle_cpu)/100;
        end
        
        function mem = get_total_mem_B()
            [~, mem_str] = system("vmstat -s -S K | grep total\ memory | awk '{print $1}'");
            mem = str2double(mem_str)*2^10;
        end
        
        function frac = get_free_mem_frac()
            frac = NN.get_free_mem_B/NN.get_total_mem_B;
        end
        
        function mem = get_free_swap_B()
            [~, mem_str] = system("vmstat -s -S K | grep free\ swap | awk '{print $1}'");
            mem = str2double(mem_str)*2^10;
        end
        
        function mem_B = get_required_mem_B(sz, dataType)
            % sz is the size vector of the matrix
            
            n = prod(sz);
            
            switch dataType
                case 'double'
                    bytesPerElement = 8;
                case 'single'
                    bytesPerElement = 4;
                case 'int32'
                    bytesPerElement = 4;
                case 'int16'
                    bytesPerElement = 2;
                case 'int8'
                    bytesPerElement = 1;
                otherwise
                    error('Unsupported data type');
            end
            
            % Calculate the total memory required in bytes
            mem_B = n*bytesPerElement;
        end
        
        function [frac_used, free_mem_B] = get_nvidia_gpu_used_memory_frac()
            [~, gpu_str] = system("nvidia-smi | sed -n '10p'");
            gpu_str = strsplit(gpu_str);
            used_mem_MB = gpu_str{9};
            MiB_loc = strfind(used_mem_MB, "MiB");
            if isempty(MiB_loc)
                error("Check processing of output of nvidia-smi")
            end
            used_mem_B = str2double(used_mem_MB(1:MiB_loc - 1))*2^20;
            total_mem_MB = gpu_str{11};
            MiB_loc = strfind(total_mem_MB, "MiB");
            if isempty(MiB_loc)
                error("Check processing of output of nvidia-smi")
            end
            total_mem_B = str2double(total_mem_MB(1:MiB_loc - 1))*2^20;
            free_mem_B = total_mem_B - used_mem_B;
            frac_used = used_mem_B/total_mem_B;
        end
        
    end
    
    
end

