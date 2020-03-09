classdef FFNNS < handle
    % FFNNS Class is a new feedforward network class used to replace the old
    % FFNN in the future
    % reachability analysis methods: 'exact-star', 'exact-polyhedron',
    % 'approx-star', 'approx-zono' (zonotope), 'abs-dom' (abstract domain)
    % Author: Dung Tran
    % Date: 27/2/2019
    
    properties
        Name = 'net';
        Layers = []; % An array of Layers, eg, Layers = [L1 L2 ...Ln]
        nL = 0; % number of Layers
        nN = 0; % number of Neurons
        nI = 0; % number of Inputs
        nO = 0; % number of Outputs
        
        % properties for reach set computation
        
        reachMethod = 'exact-star';    % reachable set computation scheme, default - 'star'
        reachOption = []; % parallel option, default - non-parallel computing
        numCores = 1; % number of cores (workers) using in computation
        inputSet = [];  % input set
        reachSet = [];  % reachable set for each layers
        outputSet = []; % output reach set
        reachTime = []; % computation time for each layers
        numReachSet = []; % number of reach sets over layers
        totalReachTime = 0; % total computation time       
        numSamples = 2000; % default number of samples using to falsify a property
        unsafeRegion = []; % unsafe region of the network
        
        
        Operations = []; % flatten a network into a sequence of operations
        
    end
    
    
    methods % constructor, evaluation, sampling, print methods
        
        % constructor
        function obj = FFNNS(varargin)
            
            switch nargin
                case 1
                    Layers = varargin{1};
                    name = 'net';
                case 2 
                    Layers = varargin{1};
                    name = varargin{2};
                otherwise
                    error('Invalid number of inputs');
            end
            
            nL = size(Layers, 2); % number of Layer
            for i=1:nL
                L = Layers(i);
                if ~isa(L, 'LayerS')
                    error('Element %d of Layers array is not a LayerS object', i);
                end
            end
            
            % check consistency between layers
            for i=1:nL-1
                if size(Layers(i).W, 1) ~= size(Layers(i + 1).W, 2)
                    error('Inconsistent dimensions between Layer %d and Layer %d', i, i + 1);
                end
            end
            
            obj.Layers = Layers;
            obj.nL = nL;    % number of layers
            obj.nI = size(Layers(1).W, 2); % number of inputs
            obj.nO = size(Layers(nL).W, 1); % number of outputs
            obj.Name = name;
            
            nN = 0;
            for i=1:nL
                nN = nN + Layers(i).N;
            end
            obj.nN = nN; % number of neurons
            
        end
        
        
        % Evaluation of a FFNN
        function y = evaluate(obj, x)
            % Evaluation of this FFNN
            % @x: input vector x
            % @y: output vector y
            
            y = x; 
            for i=1:obj.nL
                y = obj.Layers(i).evaluate(y);
            end
        
        end
        
        % Sample of FFNN
        function Y = sample(obj, V)
            % sample the output of each layer in the FFNN based on the
            % vertices of input set I, this is useful for testing.
            % @V : array of vertices to evaluate
            % @Y : output which is a cell array
        
            In = V;
            for i=1:obj.nL
                In = obj.Layers(i).sample(In);
            end          
            Y = In;
        end
        
        % check if all activation functions are piece-wise linear
        function b = isPieceWiseNetwork(obj)
            % author: Dung Tran
            % date: 9/30/2019
            
            n = obj.nL; 
            
            b = 1; 
            for i=1:obj.nL
                f = obj.Layers(i).f;                
                if ~strcmp(f, 'poslin') && ~strcmp(f, 'purelin') && ~strcmp(f, 'satlin') && ~strcmp(f, 'satlins')
                    b = 0;
                    return;
                end
                              
            end
            
        end
        
        
        % print information to a file
        function print(obj, file_name)
            % @file_name: name of file you want to store all data
            % information
            f = fopen(file_name, 'w');
            fprintf(f, 'Feedforward Neural Network Information\n');
            fprintf(f, '\nNumber of layers: %d', obj.nL);
            fprintf(f, '\nNumber of neurons: %d', obj.nN);
            fprintf(f, '\nNumber of inputs: %d', obj.nI);
            fprintf(f, '\nNumber of outputs: %d', obj.nO);
            
            if ~isempty(obj.reachSet)
                fprintf(f, '\n\nReach Set Information');
                fprintf(f, '\nReachability method: %s', obj.reachMethod);
                fprintf(f, '\nNumber of cores used in computation: %d', obj.numCores);
                
                for i=1:length(obj.reachSet)-1
                    fprintf(f, '\nLayer %d reach set consists of %d sets that are computed in %.5f seconds', i, obj.numReachSet(i), obj.reachTime(i));
                end
                fprintf(f, '\nOutput Layer reach set consists of %d sets that are computed in %.5f seconds', obj.numReachSet(obj.nL), obj.reachTime(obj.nL)); 
                fprintf(f, '\nTotal reachable set computation time: %.5f', obj.totalReachTime);
                
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
        
        function [S, t] = reach(varargin)            
            % @I: input set, a star set          
            % @method: = 'exact-star' or 'approx-star' -> compute reach set using stars
            %            'abs-dom' -> compute reach set using abstract
            %            domain (support in the future)
            %            'face-latice' -> compute reach set using
            %            face-latice (support in the future)
            
            % @R: output set 
            % @t : computation time 
            
            % @numOfCores: number of cores you want to run the reachable
            % set computation, @numOfCores >= 1, maximum is the number of
            % cores in your computer. We also can run this method on a
            % clusters. You need to change: parpool('local', numOfCores) to
            % parcluster(your_cluster_profile), you also need to set up
            % your local clusters with installed MPT toolbox also. see: 
            % https://www.mathworks.com/help/distcomp/discover-clusters-and-use-cluster-profiles.html
           
            
            % parse inputs 
            switch nargin
                case 4
                    obj = varargin{1}; % FFNNS object
                    obj.inputSet = varargin{2}; % input set
                    obj.reachMethod = varargin{3}; % reachability analysis method
                    obj.numCores = varargin{4}; % number of cores used in computation
                    
                case 3
                    obj = varargin{1};
                    obj.inputSet = varargin{2};
                    obj.reachMethod = varargin{3};
                    obj.numCores = 1;
                    
                case 2
                    obj = varargin{1};
                    obj.inputSet = varargin{2};
                    obj.reachMethod = 'exact-star';
                    obj.numCores = 1;
                    
                otherwise
                    error('Invalid number of input arguments (should be 1, 2 or 3)');
            end
            
            
            % if reachability analysis method is an over-approximate
            % method, we use 1 core for computation
            if ~strcmp(obj.reachMethod, 'exact-star') && ~strcmp(obj.reachMethod, 'exact-polyhedron')
                obj.numCores = 1;
            end
            
            % Zonotope method accepts both star and zonotope input set
            if strcmp(obj.reachMethod, 'approx-zono') && isa(varargin{2}, 'Star')
                obj.inputSet = varargin{2}.getZono;
            end
            
            if obj.numCores == 1
                obj.reachOption = []; % don't use parallel computing
            else
                obj.start_pool;       % start parallel pool in Matlab
                obj.reachOption = 'parallel';
            end
            
            obj.reachSet = cell(1, obj.nL);
            obj.numReachSet = zeros(1, obj.nL);
            obj.reachTime = [];
                        
            % compute reachable set
            In = obj.inputSet;
                        
            for i=1:obj.nL               
                
                fprintf('\nComputing reach set for Layer %d ...', i);
                
                % estimate reachability analysis time from the second layer
                if i > 1
                    nP = length(obj.reachSet{1, i - 1});
                    if i==2
                        nP1 = 1;
                    else
                        nP1 = length(obj.reachSet{1, i-2});
                    end
                    rt = obj.reachTime(i-1);
                    
                    estimatedTime = rt*(nP/nP1)*(obj.Layers(i).N / obj.Layers(i-1).N);
                    
                    fprintf('\nEstimated computation time: ~ %.5f seconds', estimatedTime);
                end
                
                st = tic;
                In = obj.Layers(i).reach(In, obj.reachMethod, obj.reachOption);
                t1 = toc(st);
                
                %obj.reachSet{1, i} = In;
                obj.numReachSet(i) = length(In);
                obj.reachTime = [obj.reachTime t1];
                
                fprintf('\nExact computation time: %.5f seconds', t1);               
                fprintf('\nNumber of reachable set at the output of layer %d: %d', i, length(In));               
                               
            end
            
            obj.outputSet = In;
            obj.totalReachTime = sum(obj.reachTime);
            S = obj.outputSet;
            t = obj.totalReachTime;                      
            fprintf('\nTotal reach set computation time: %.5f seconds', obj.totalReachTime);
            fprintf('\nTotal number of output reach sets: %d', length(obj.outputSet));
            
        end
        
        
    end
    
    methods
        
        function [safe, vt, counterExamples] = verify(varargin)            
            % 1: @I: input set, need to be a star set
            % 2: @U: unsafe region, a set of HalfSpaces
            % 3: @method: = 'star' -> compute reach set using stars
            %            'abs-dom' -> compute reach set using abstract
            %            domain (support in the future)
            %            'face-latice' -> compute reach set using
            %            face-latice (support in the future)
            % 4: @numOfCores: number of cores you want to run the reachable
            % set computation, @numOfCores >= 1, maximum is the number of
            % cores in your computer.
            
            % 5: @n_samples : number of simulations used for falsification if
            % using over-approximate reachability analysis, i.e.,
            % 'approx-zono' or 'abs-dom' or 'abs-star'
            % note: n_samples = 0 -> do not do falsification
            
            switch nargin
                case 3
                    obj = varargin{1};
                    obj.inputSet = varargin{2};
                    obj.unsafeRegion = varargin{3};
                    
                case 4
                    obj = varargin{1};
                    obj.inputSet = varargin{2};
                    obj.unsafeRegion = varargin{3};
                    obj.reachMethod = varargin{4};
                    
                case 5
                    obj = varargin{1};
                    obj.inputSet = varargin{2};
                    obj.unsafeRegion = varargin{3};
                    obj.reachMethod = varargin{4};
                    obj.numCores = varargin{5};
                case 6
                    obj = varargin{1};
                    obj.inputSet = varargin{2};
                    obj.unsafeRegion = varargin{3};
                    obj.reachMethod = varargin{4};
                    obj.numCores = varargin{5};
                    obj.numSamples = varargin{6};
                 
                otherwise
                    error('Invalid number of inputs, should be 2, 3, 4 or 5');
            end
            
            t = tic; 
            fprintf('\nPerform fasification with %d random simulations', obj.numSamples);
            counterExamples = obj.falsify(obj.inputSet, obj.unsafeRegion, obj.numSamples);
            
            if ~isempty(counterExamples)
                safe = 0;
            else
                fprintf('\nNo counter examples found, verify the safety using reachability analysis');
                % perform reachability analysis
                [R,~] = obj.reach(obj.inputSet, obj.reachMethod, obj.numCores);   
                
                if strcmp(obj.reachMethod, 'exact-star')
                    
                    n = length(R);
                    counterExamples = [];
                    
                    for i=1:n
                        if ~isempty(R(i).intersectHalfSpace(obj.unsafeRegion.G, obj.unsafeRegion.g))
                            counterExamples = [counterExamples Star(obj.inputSet.V, R(i).C, R(i).d, R(i).predicate_lb, R(i).predicate_ub)];
                        end
                    end
                    
                    if isempty(counterExamples)
                        safe = 1;
                    else
                        safe = 0;
                    end
                    
                    
                else
                    if strcmp(obj.reachMethod, 'zono')
                        R = R.toStar;
                    end
                    
                    if isempty(R.intersectHalfSpace(obj.unsafeRegion.G, obj.unsafeRegion.g))
                        safe = 1;
                    else
                        safe = 2;
                    end
                    
                end
                
            end
            
            vt = toc(t);
         
        end
    end
    
    methods % Input to Output Sensitivity
        
        function [maxSensInputId, maxSensVal] = getMaxSensitiveInput(obj, I, output_mat)
            % @I: input, is a star set, or box, or zono
            % @output_mat: output matrix, y = C*x, x is the output vector
            % of the networks
            % @maxSensInputId: the id of the input that is most sensitive
            % with the output changes 
            % @maxSensVal: percentage of sensitivity of all inputs
            
            % author: Dung Tran
            % date: 10/19/2019
            
            
            if ~isa(I, 'Box')
                I1 = I.getBox;
            else
                I1 = I;
            end
            
            if I1.dim ~= obj.nI
                error('Inconsistency between the input set and the number of inputs in the networks');
            end
            
            if obj.nO ~= size(output_mat, 2)
                error('Inconsistency between the output matrix and the number of outputs of the networks');
            end
            
            maxSensVal = zeros(1, obj.nI);
            [R, ~] = obj.reach(I1.toZono, 'approx-zono');
            R = R.affineMap(output_mat, []);
            B = R.getBox;
            max_rad = max(B.ub);
            
            for i=1:obj.nI
                I2 = I1.singlePartition(i, 2); % partition into two boxes at index i
                [R2, ~] = obj.reach(I2(1).toZono, 'approx-zono');
                R2 = R2.affineMap(output_mat, []);
                B2 = R2.getBox;
                max_rad2 = max(B2.ub); 
                maxSensVal(i) = (max_rad - max_rad2)/max_rad;
            end
            
            [~,maxSensInputId] = max(maxSensVal);
   
        end
        
        
        function Is = partitionInput_MSG(obj, I, U)
            % @I: input, is a star set, or box, or zono
            % @U: unsafe region, a HalfSpace object
            % @Is: two partitioned inputs, the last one is close to the
            %      unsafe region
            
            % author: Dung Tran
            % date: 10/19/2019
            
            
            if ~isa(I, 'Box')
                I1 = I.getBox;
            else
                I1 = I;
            end
            
            if I1.dim ~= obj.nI
                error('Inconsistency between the input set and the number of inputs in the networks');
            end
            
            if ~isa(U, 'HalfSpace')
                error('Unsafe region is not a HalfSpace object');
            end
            
            if obj.nO ~= size(U.G, 2)
                error('Inconsistency between the output matrix and the number of outputs of the networks');
            end
            
            maxSensVal = zeros(1, obj.nI);
            [R, ~] = obj.reach(I1.toZono, 'approx-zono');
            R = R.affineMap(U.G, -U.g);
            B = R.getBox;
            max_rad = max(B.ub);
            
            for i=1:obj.nI
                I2 = I1.singlePartition(i, 2); % partition into two boxes at index i
                [R2, ~] = obj.reach(I2(1).toZono, 'approx-zono');
                R2 = R2.affineMap(U.G, -U.g);
                B2 = R2.getBox;
                max_rad2 = max(B2.ub); 
                maxSensVal(i) = (max_rad - max_rad2)/max_rad;
            end
            
            [~,maxSensInputId] = max(maxSensVal);
            
            I1 = I.partition(maxSensInputId, 2); 
            R1 = obj.reach(I1(1).toZono, 'approx-zono');
            R1 = R1.affineMap(U.G, -U.g);
            B1 = R1.getBox;
            max_y1 = max(B1.ub);
            R2 = obj.reach(I1(2).toZono, 'approx-zono');
            R2 = R2.affineMap(U.G, -U.g);
            B2 = R2.getBox;
            max_y2 = max(B2.ub);
            
            if max_y1 <= max_y2
                Is = [I1(2) I1(1)];
            else
                Is = [I1(1) I1(2)];
            end
            
            
            
   
        end
        
        
        % depth first search for max sensitive inputs in partitioning
        function [maxSensInputs, maxSensVals] = searchMaxSensInputs(obj, I, output_mat, k, maxSen_lb)
            % @I: input, is a star set, or box, or zono
            % @output_mat: output matrix, y = C*x, x is the output vector
            %              of the networks
            % @k: depth of search tree
            % @maxSen_lb: the search is stop if the max sensitive value
            %             is smaller than the maxSen_lb
            
            % @maxSensInputs: a sequence of input indexes that are most sensitive
            %                 with the output changes 
            % @maxSensVals: percentage of sensitivity corresponding to
            %               above sequence of inputs
          
            
            
            % author: Dung Tran
            % date: 10/19/2019
            
            if k < 1 
                error('Invalid depth of search tree');
            end
            
            if ~isa(I, 'Box')
                error('Input is not a box');
            end
            
            maxSensInputs = [];
            maxSensVals = [];
            for i=1:k
                
                if i==1
                    
                    [max_id, sens_vals] = obj.getMaxSensitiveInput(I, output_mat);
                    maxSensInputs = [maxSensInputs max_id];
                    maxSensVals = [maxSensVals max(sens_vals)];
                    I1 = I; 
                else
                    I2 = I1.singlePartition(maxSensInputs(i-1), 2);
                    I2 = I2(1); 
                    [max_id, sens_vals] = obj.getMaxSensitiveInput(I2, output_mat);
                    maxSensInputs = [maxSensInputs max_id];
                    maxSensVals = [maxSensVals max(sens_vals)];
                    I1 = I2;
                end
                
                if ~isempty(maxSen_lb) && (maxSensVals(i) < maxSen_lb)
                    maxSensInputs(i) = [];
                    maxSensVals(i) = [];
                    break;
                end

            end
              
        end
        
        
        % verify safety property using max sensitivity guided (MSG) method       
        function [safe, VT, counterExamples] = verify_MSG(obj, I, reachMethod, k, sens_lb, U)
            % @I: input set, a box
            % @reachMethod: reachability method
            % @k: depth of search tree for max sensitive inputs
            % @sens_lb: lowerbound of sensitive value
            % @U: unsafe region, a halfspace object
            % @safe: = 1: safe
            %        = 2: unknown
            %        = 0: unsafe
           
            % author: Dung Tran
            % date: 10/20/2019
                       
            t = tic; 
            [maxSensInputs, ~] = obj.searchMaxSensInputs(I, U.G, k, sens_lb);
            n = length(maxSensInputs);
            Is = I.partition(maxSensInputs, 2*ones(1, n));           
            I1 = []; % set of partitioned inputs
            N = 2^n; % number of partitioned inputs
            for i=1:N                
                if strcmp(reachMethod, 'approx-zono')
                    I1 = [I1 Is(i).toZono];
                elseif strcmp(reachMethod, 'approx-star') || strcmp(reachMethod, 'abs-dom')
                    I1 = [I1 Is(i).toStar];
                else
                    error('reachmethod should be approx-zono or approx-star or abs-dom');
                end
            end
            
            safe_vec = zeros(1, N);
            [R, ~] = obj.reach(I1, reachMethod); % perform reachability anlaysis
            for i=1:N
                S = R(i).intersectHalfSpace(U.G, U.g);
                if isempty(S)
                    safe_vec(i) = 1;
                else
                    counterExamples = obj.falsify(I1(i), U, 1);
                    if length(counterExamples) >= 1
                        safe_vec(i) = 0;
                        break;
                    else
                        safe_vec(i) = 2;
                    end
                end
            end
            
            if sum(safe_vec) == N
                safe = 'SAFE';
                counterExamples = [];
                fprintf('\nTHE NETWORK IS SAFE');
            else
                
                if ~isempty(counterExamples)
                    safe = 'UNSAFE';
                    fprintf('\nTHE NETWORK IS UNSAFE');
                else
                    safe = 'UNKNOWN';
                    fprintf('\nTHE SAFETY OF THE NETWORK IS UNKNOWN');
                end

            end
            
            VT = toc(t);
            
        end
        
        
        % get counter input candidate after a single partition using MSG
        function counterInputCand = getStepCounterInputCand_MSG(obj, I, U)
            % @I: input, is a star set, or box, or zono
            % @output_mat: output matrix, y = C*x, x is the output vector
            %              of the networks
            % @k: 
            % @maxSensInputId: the id of the input that is most sensitive
            %                  with the output changes
            % @counterInputCand: counter input candidate
            
            % author: Dung Tran
            % date: 10/20/2019
            
            
            if ~isa(I, 'Box')
                I1 = I.getBox;
            else
                I1 = I;
            end
            
            if I1.dim ~= obj.nI
                error('Inconsistency between the input set and the number of inputs in the networks');
            end
            
            if obj.nO ~= size(U.G, 2)
                error('Inconsistency between the unsafe region and the number of outputs of the networks');
            end
            
            maxSensVal = zeros(1, obj.nI);

            [R, ~] = obj.reach(I1.toZono, 'approx-zono');
            R = R.affineMap(U.G, -U.g);
            B = R.getBox;
            max_rad = max(B.ub - B.lb);
            
            for i=1:obj.nI
                I2 = I1.singlePartition(i, 2); % partition into two boxes at index i
                [R1, ~] = obj.reach(I2(1).toZono, 'approx-zono');
                R1 = R1.affineMap(U.G, -U.g);
                B1 = R1.getBox;                             
                max_rad1 = max(B1.ub - B1.lb); 
                maxSensVal(i) = (max_rad - max_rad1)/max_rad;
            end
            
            [~,maxSensInputId] = max(maxSensVal);
            
            I2 = I1.singlePartition(maxSensInputId, 2);
            [R1, ~] = obj.reach(I2(1).toZono, 'approx-zono');
            R1 = R1.affineMap(U.G, -U.g);
            B1 = R1.getBox;                             
            max_y1 = max(B1.ub);
            
            [R2, ~] = obj.reach(I2(2).toZono, 'approx-zono');
            R2 = R2.affineMap(U.G, -U.g);
            B2 = R2.getBox;                             
            max_y2 = max(B2.ub);
            
            if max_y1 <= max_y2
                counterInputCand = I2(1);
            else
                counterInputCand = I2(2);
            end

            
        end
        
        
        % Depth First Search for Counter Input Candidate
        function counterInputCand = searchCounterInputCand_MSG(obj, I, U, k)
            % @I: input, is a star set, or box, or zono
            % @U: unsafe region of the networks
            % @k: depth of search tree
            % @maxSensInputId: the id of the input that is most sensitive
            %                  with the output changes
            % @counterInputCand: counter input candidate
            
            % author: Dung Tran
            % date: 10/20/2019
            
            if k < 1
                error('depth of search tree should be >= 1');
            end
            
            counterInputCand = I; 
            for i=1:k
                counterInputCand = obj.getStepCounterInputCand_MSG(counterInputCand, U);
            end

        end
        
        
        % Depth First Seach for Falsification using Maximum Sensitive
        % Guided Method
        
        function  counterInputs = falsify_MSG(obj, I, U)
            % @I: input set, is a box          
            % @U: unsafe region, a halfspace object
            % @counterExamples: counter example inputs
            
            % author: Dung Tran
            % date: 10/20/2019
            
            
            
        end
        
        
        function [safe, VT, counterInputs] = verify_MSG2(obj, I, U)
            % @I: input set, is a box          
            % @U: unsafe region, a halfspace object
            % @counterExamples: counter example inputs
            
            % author: Dung Tran
            % date: 10/20/2019
            
            t = tic;
            if ~isa(I, 'Box')
                error('Input is not a Box');
            end
            
            if I.dim ~= obj.nI
                error('Inconsistency between input set and the number of network inputs');
            end
            
            if ~isa(U, 'HalfSpace')
                error('Unsafe region is not a HalfSpace Class');
            end
            
            if size(U.G, 2) ~= obj.nO
                error('Inconsistent between the Unsafe region and the number of outputs of the network');
            end
            
            I0 = I; % a queue to store number of partitioned input sets
            safe = 'SAFE';
            counterInputs = [];
            while ~isempty(I0) 
                n = length(I0); 
                I1 = I0(n); % pop the last input set in the queue
                I0(n) = []; 
                [R, ~] = obj.reach(I1.toZono, 'approx-zono');
                S = R.intersectHalfSpace(U.G, U.g);                
                if ~isempty(S)                    
                    y1 = obj.evaluate(I1.lb);
                    y2 = obj.evaluate(I1.ub);
                    y3 = obj.evaluate(0.5*(I1.lb + I1.ub));                  
                    ys = [];
                    if U.contains(y1)
                        ys = [ys y1];
                    end
                    if U.contains(y2)
                        ys = [ys y2];
                    end
                    if U.contains(y3)
                        ys = [ys y3];
                    end                    
                    if ~isempty(ys)
                        safe = 'UNSAFE';
                        counterInputs = ys;
                        break;
                    else
                        I0 = [I0 obj.partitionInput_MSG(I1, U)]; % put new partitioned input into the queue
                    end 
                end
                
                

            end
            
            VT = toc(t);
            
            
        end
        
    end
    
    
    
    
    methods % checking safety method or falsify safety property
        
        function [safe, t, counter_inputs] = isSafe(varargin)
            % 1: @I: input set, need to be a star set
            % 2: @U: unsafe region, a set of HalfSpaces
            % 3: @method: = 'star' -> compute reach set using stars
            %            'abs-dom' -> compute reach set using abstract
            %            domain (support in the future)
            %            'face-latice' -> compute reach set using
            %            face-latice (support in the future)
            % 4: @numOfCores: number of cores you want to run the reachable
            % set computation, @numOfCores >= 1, maximum is the number of
            % cores in your computer.
            
            % 5: @n_samples : number of simulations used for falsification if
            % using over-approximate reachability analysis, i.e.,
            % 'approx-zono' or 'abs-dom' or 'abs-star'
            % note: n_samples = 0 -> do not do falsification
            
                       
            % @safe: = 1 -> safe, = 0 -> unsafe, = 2 -> uncertain 
            % @t : verification time 
            % @counter_inputs
            
            % author: Dung Tran
            % date: 3/27/2019
            
            % parse inputs 
            switch nargin
                
                case 6
                    
                    obj = varargin{1}; % FFNNS object
                    I = varargin{2};
                    U = varargin{3};
                    method = varargin{4};
                    n_samples = varargin{5}; 
                    numOfCores = varargin{6};
                    
                    
                 case 5
                    
                    obj = varargin{1}; % FFNNS object
                    I = varargin{2};
                    U = varargin{3};
                    method = varargin{4};
                    n_samples = varargin{5}; 
                    numOfCores = 1;
                    
                                                    
                case 4
                    
                    obj = varargin{1}; % FFNNS object
                    I = varargin{2};
                    U = varargin{3};
                    method = varargin{4};
                    n_samples = 1000; % using 1000 samples
                    numOfCores = 1;
     
                case 3
                    
                    obj = varargin{1}; % FFNNS object
                    I = varargin{2};
                    U = varargin{3};
                    method = 'exact-star';
                    n_samples = 0; 
                    numOfCores = 1;
                    
                otherwise
                    error('Invalid number of input arguments (should be 3, 4, 5 or 7), please look at ../tests/nn/ffnns/ for usages of this method');
            end
                              
            start = tic;
                        
            if isempty(U)
                error('Please specify unsafe region using Half-space class');
            end
            
            % performing reachability analysis
            [R, ~] = obj.reach(I, method, numOfCores);        
                      
            % check safety
            n = length(R);
            m = length(U); 
            R1 = [];
            for i=1:n
                if isa(R(i), 'Zono')
                    B = R(i).getBox;
                    R1 = [R1 B.toStar]; % transform to star sets;
                else
                    R1 = [R1 R(i)];
                end
            end
                                    
            violate_inputs = [];            
            
            if numOfCores == 1

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

            elseif numOfCores > 1

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
            
            
            if isempty(violate_inputs)
                
                safe = 1;  
                counter_inputs = []; 
                fprintf('\nThe network is safe');
                
            else

                if strcmp(method, 'exact-star')
                    
                    safe = 0;  
                    counter_inputs = violate_inputs; % exact-method return complete counter input set
                    fprintf('\nThe network is unsafe, counter inputs contains %d stars', length(counter_inputs));
                    
                else
                    
                    if n_samples == 0
                        fprintf('\nDo not do falsification since n_samples = 0, you can choose to do falsification by set n_samples value > 0');
                        safe = 2;
                        counter_inputs = [];
                    else
                        
                        counter_inputs = obj.falsify(I, U, n_samples);
                    
                        if isempty(counter_inputs)
                            safe = 2;
                            fprintf('\nSafety is uncertain under using %d samples to falsify the network', n_samples);
                            fprintf('\nYou can try to increase the samples for finding counter inputs');
                        else
                            safe = 0;
                            fprintf('\nThe network is unsafe, %d counter inputs are found using %d simulations', length(counter_inputs), n_samples);
                        end

                        
                    end
                    
                    
                end


            end

            t = toc(start);
             
            
        end
                
        
        % falsify safety property using random simulation        
        function counter_inputs = falsify(obj, I, U, n_samples)
            % @input: star set input
            % @U: unsafe region, a set of HalfSpaces
            % @n_samples: number of samples used in falsification
            % @counter_inputs: counter inputs that falsify the property
            
            % author: Dung Tran
            % date: 3/27/2019
            
            counter_inputs = [];
            
            if isa(I, 'Zono') || isa(I, 'Box')
                I1 = I.toStar;
            elseif isa(I, 'Star')
                I1 = I;
            else
                error('Unknown set representation');
            end
                        
            m = length(U);
            
            for i=1:m
                if ~isa(U(i), 'HalfSpace')
                    error('%d^th unsafe region is not a HalfSpace', i);
                end
            
            end
            
            if n_samples < 1
                error('Invalid number of samples');
            end
            
            V = I1.sample(n_samples);
            
            n = size(V, 2); % number of samples 
           
            for i=1:n
                
                y = obj.evaluate(V(:, i));
                
                for j=1:m
                    if U(j).contains(y)
                        counter_inputs = [counter_inputs V(:, i)];
                    end
                end
                
            end
            
            
        end
        
         
    end
    
    
    methods % checking robustness and get robustness bound of feedforward networks
        
        
        % Problem description:
        % checking robustness of FFNN corresponding to an input and
        % L_infinity norm bound disturbance
        
        % x is input, x' is a disturbed input such that ||x'- x|| <= e
        % y = FFNN(x), and y' = FFNN(x')
        % FFNN is called robust if y' is in a robust region RB(y') defined by user
        % for example, let F is image classification network
        % inputs of F is an 28x28 = 784 pixel images, F has one output with
        % range between 0 to 9 to reconize 10 digits, 0, 1, 2 ..., 9
        
        % let consider digit one as an input x, we disturb digit one with
        % some bounded disturabnce e, we have ||x' - x|| <= e
        
        % to be called robust, F need produce the output sastify 0.5 <= y' <= 1.5
        % If output is violated this specification, F may misclasify digit
        % one. Thefore the un-robust regions in this case is S1 = y' < 0.5 and S2 = y'> 1.5
        
        % If robustness is not guaranteed, the method search for some
        % adverserial examples, i.e., find some x' that makes F violate the
        % robustness property.
        
        function [robust, t, adv_inputs] = isRobust(varargin)
            % 1: @input_vec: input vector x
            % 2: @dis_bound: disturbance bound
            % 3: @un_robust_reg: un-robust-region, a star
            % 4: @method: = 'exact-star' or 'approx-star' or 'approx-zono' or 'abs-dom'
            % 5: @lb_allowable: allowable lower bound of disturbed input:    lb_allowable(i) <= x'[i]
            % 6: @ub_allowable: allowable upper bound of disturbed output:    ub_allowable(i) >= x'[i] 
            % x' is the disturbed vector by disturbance bound, |x' - x| <= dis_bound
            % x'[i] >= lb_allowable, x'[i] <= ub_allowable[i]
            % 7: @n_samples: number of samples used to find counter examples
            % if using over-approximate reachability analysis methods
            % 8: @numCores:  number of cores used in computation
            
            
            % @robust: = 1-> robust
            %        : = 0 -> unrobust
            %        : = 2 -> uncertain, cannot find counter example
            % @adv_inputs: adverserial inputs
            
            % author: Dung Tran
            % date: 3/27/2019
            
            
            start = tic; 
            % parse inputs
            switch nargin
                case 9
                    obj = varargin{1}; % FFNNS object
                    input_vec = varargin{2}; % input vec
                    dis_bound = varargin{3}; % disturbance bound
                    un_robust_reg = varargin{4}; % un-robust region
                    method = varargin{5}; % reachability analysis method
                    lb_allowable = varargin{6}; % allowable lower bound on disturbed inputs
                    ub_allowable = varargin{7}; % allowable upper bound on disturbed inputs
                    n_samples = varargin{8}; % number of samples used for finding counter examples
                    num_cores = varargin{9}; % number of cores used in computation
                    
                    % check consistency
                    if ~isempty(lb_allowable) && ((length(lb_allowable) ~= length(ub_allowable)) || (length(lb_allowable) ~= length(input_vec)))
                        error('Inconsistent dimensions between allowable lower-, upper- bound vectors and input vector');
                    end
                    
                case 8
                    
                    obj = varargin{1}; % FFNNS object
                    input_vec = varargin{2}; % input vec
                    dis_bound = varargin{3}; % disturbance bound
                    un_robust_reg = varargin{4}; % un-robust region
                    method = varargin{5}; % reachability analysis method
                    lb_allowable = varargin{6}; % allowable lower bound on disturbed inputs
                    ub_allowable = varargin{7}; % allowable upper bound on disturbed inputs
                    n_samples = varargin{8}; % number of samples used for finding counter examples
                    num_cores = 1; % number of cores used in computation
                    
                case 6
                    
                    obj = varargin{1}; % FFNNS object
                    input_vec = varargin{2}; % input vec
                    dis_bound = varargin{3}; % disturbance bound
                    un_robust_reg = varargin{4}; % un-robust region
                    method = varargin{5}; % reachability analysis method
                    lb_allowable =[]; % allowable lower bound on disturbed inputs
                    ub_allowable = []; % allowable upper bound on disturbed inputs
                    n_samples = varargin{6}; % number of samples used for finding counter examples
                    num_cores = 1; % number of cores used in computation
                    
                case 5
                    
                    obj = varargin{1}; % FFNNS object
                    input_vec = varargin{2}; % input vec
                    dis_bound = varargin{3}; % disturbance bound
                    un_robust_reg = varargin{4}; % un-robust region
                    method = varargin{5}; % reachability analysis method
                    lb_allowable =[]; % allowable lower bound on disturbed inputs
                    ub_allowable = []; % allowable upper bound on disturbed inputs
                    n_samples = 1000; % use 1000 number of samples used for finding counter examples
                    num_cores = 1; % number of cores used in computation
                    
                    
                otherwise
                    error('Invalid number of input arguments (should be 4 or 5 or 8), please look at ../tests/nn/ffnns/ for usages of this method');
            end
            
            
            % construct input set
            n = length(input_vec);
            lb = input_vec;
            ub = input_vec;
            
            if nargin == 8 || nargin == 9
                
                for i=1:n
                    
                    if lb(i) - dis_bound > lb_allowable(i)
                        lb(i) = lb(i) - dis_bound;
                    else
                        lb(i) = lb_allowable(i);
                    end
                    if ub(i) + dis_bound < ub_allowable(i)
                        ub(i) = ub(i) + dis_bound;
                    else
                        ub(i) = ub_allowable(i);
                    end  
                    
                end
                
            else    
                for i=1:n    
                    lb(i) = lb(i) - dis_bound;
                    ub(i) = ub(i) + dis_bound;
                end
            end
            
            
            % input set to check robustness
            I = Star(lb, ub);
                                   
            [robust, ~, adv_inputs] = obj.isSafe(I, un_robust_reg, method, n_samples, num_cores);
            
            if robust == 1
                fprintf('\nThe network is robust with the disturbance dis_bound = %.5f', dis_bound);
            elseif robust == 0
                fprintf('\nThe network is not robust with the disturbance dis_bound = %.5f, counter examples are found', dis_bound);
            elseif robust == 2
                fprintf('\nThe robustness of the network is uncertain with the disturbance dis_bound = %.5f, we cannot find counter examples, you can try again with a larger n_samples', dis_bound);
            end
                            
            t = toc(start);
                       
        end
        
        
        % find maximum robustness value, i.e., maximum disturbance bound
        % that the network is still robust
        
        function [robustness_bound, t] = get_robustness_bound(varargin)
            % 1: @input_vec: input point
            % 2: @init_dis_bound: initial disturbance bound
            % 3: @dis_bound_step: a step to increase/decrease disturbance bound
            % 4: @max_steps: maximum number of steps for searching
            % 5: @lb_allowable: allowable lower bound of disturbed input:    lb_allowable(i) <= x'[i]
            % 6: @ub_allowable: allowable upper bound of disturbed output:    ub_allowable(i) >= x'[i] 
            % x' is the disturbed vector by disturbance bound, |x' - x| <= dis_bound
            % x'[i] >= lb_allowable, x'[i] <= ub_allowable[i]
            % 7: @method: = 'exact-star' or 'approx-star' or 'approx-zono' or 'abs-dom'
            % 8: @n_samples: number of samples used to find counter examples
            % if using over-approximate reachability analysis methods
            % 9: @numCores:  number of cores used in computation
            
            % @robustness_bound: robustness bound w.r.t @input_vec
            % @t: computation time
            
            % author: Dung Tran
            % date: 3/29/2019
            
            start = tic;
            switch nargin
                
                case 11
                    obj = varargin{1}; % FFNNS object
                    input_vec = varargin{2}; % input vec
                    init_dis_bound = varargin{3}; % initial disturbance bound for searching
                    tolerance = varargin{4}; % tolerance (accuracy) for searching
                    max_steps = varargin{5}; % maximum searching steps
                    lb_allowable = varargin{6}; % allowable lower bound on disturbed inputs
                    ub_allowable = varargin{7}; % allowable upper bound on disturbed inputs
                    un_robust_reg = varargin{8}; % un-robust region
                    method = varargin{9}; % reachability analysis method
                    n_samples = varargin{10}; % number of samples used for finding counter examples
                    num_cores = varargin{11}; % number of cores used in computation
                    
                case 9
                    
                    obj = varargin{1}; % FFNNS object
                    input_vec = varargin{2}; % input vec
                    init_dis_bound = varargin{3}; % initial disturbance bound for searching
                    tolerance = varargin{4}; % tolerance (accuracy) for searching
                    max_steps = varargin{5}; % maximum searching steps
                    lb_allowable = varargin{6}; % allowable lower bound on disturbed inputs
                    ub_allowable = varargin{7}; % allowable upper bound on disturbed inputs
                    un_robust_reg = varargin{8}; % un-robust region
                    method = varargin{9}; % reachability analysis method
                    n_samples = 1000; % number of samples used for finding counter examples
                    num_cores = 1; % number of cores used in computation
                    
                case 7
                                        
                    obj = varargin{1}; % FFNNS object
                    input_vec = varargin{2}; % input vec
                    init_dis_bound = varargin{3}; % initial disturbance bound for searching
                    tolerance = varargin{4}; % tolerance (accuracy) for searching
                    max_steps = varargin{5}; % maximum searching steps
                    lb_allowable = []; % allowable lower bound on disturbed inputs
                    ub_allowable = []; % allowable upper bound on disturbed inputs
                    un_robust_reg = varargin{6}; % un-robust region
                    method = varargin{7}; % reachability analysis method 
                    n_samples = 1000; % number of samples used for finding counter examples
                    num_cores = 1; % number of cores used in computation
                    
                case 6
                                        
                    obj = varargin{1}; % FFNNS object
                    input_vec = varargin{2}; % input vec
                    init_dis_bound = varargin{3}; % initial disturbance bound for searching
                    tolerance = varargin{4}; % tolerance (accuracy) for searching
                    max_steps = varargin{5}; % maximum searching steps
                    lb_allowable = []; % allowable lower bound on disturbed inputs
                    ub_allowable = []; % allowable upper bound on disturbed inputs
                    un_robust_reg = varargin{6}; % un-robust region
                    method = 'exact-star'; % used exact-star reachability analysis method
                    n_samples = 1000; % number of samples used for finding counter examples
                    num_cores = 1; % number of cores used in computation
                    
                    
                otherwise
                    error('Invalid number of input arguments (should be 5, 6, 8 or 10), please look at ../tests/nn/ffnns/ for usages of this method');
            end
            
            
            k = 1;
            b = init_dis_bound;
            bmax = 0;          
            while (k < max_steps)
                
                fprintf('\nSearching maximum robustness value at step k = %d...', k)
                fprintf('\nCheck robustness with disturbance bound dis_bound = %.5f...', b);
                [robust, ~, ~] = obj.isRobust(input_vec, b, un_robust_reg, method, lb_allowable, ub_allowable, n_samples, num_cores);
                
                if robust == 1
                    bmax = b;
                    b = b + tolerance;
                else
                    b = b - tolerance;
                    if b == bmax
                        break;
                    end
                end
                
                k = k + 1;
                
                
            end
            
            if k == max_steps
                fprintf('\nCannot find robustness value, increase number of searching steps, i.e., max_steps, and try again');
                robustness_bound = [];
            else
                fprintf('\nMaximum robustness value = %.5f is found at k = %d with error tolerance tol = %.5f', bmax, k, tolerance);
                robustness_bound = bmax;
            end
            
            t = toc(start);
            
        end
        
        
               
    end
    
    
    methods % verify using Deep First Search + exact star for verification (on testing phase)
        
        
        % flatten a FFNN into a sequence of operations for reachability
        function flatten(obj, reachMethod)
            % @reachMethod: reachability method
            
            % author: Dung Tran
            % date: 1/18/2020
            
            
            Ops = [];
            
            for i=1:obj.nL
                Op = obj.Layers(i).flatten(reachMethod);
                Ops = [Ops Op];
            end
            
            obj.Operations = Ops;
                       
        end
        
        
        function [safe, CEx] = verify_DFS(varargin)
            % @inputSets: a star set
            % @unsafeRegion: a HalfSpace object
            % @numCores: number of cores used for verification
            % @safe:  = 'safe' or 'unsafe' or 'unknown'
            % @CEx: counter examples
            
            % author: Dung Tran
            % date: 1/18/2020
            % update: 3/8/2020
            
            
            % passing inputs
            
            if mod(nargin, 2) == 0
                error('Invalid number of arguments');
            end
            
            obj = varargin{1};
            
            for i=2:nargin-1
                
                if mod(i, 2) == 0
                   
                    if strcmp(varargin{i}, 'InputSet')
                        obj.inputSet = varargin{i+1};                      
                    elseif strcmp(varargin{i}, 'UnsafeRegion')
                        obj.unsafeRegion = varargin{i+1};
                    elseif strcmp(varargin{i}, 'ReachMethod')
                        obj.reachMethod = varargin{i+1};
                    elseif strcmp(varargin{i}, 'NumCores')
                        obj.numCores = varargin{i+1};
                    end
                    
                end
                
            end         
            
            obj.flatten(obj.reachMethod);
            
            N = length(obj.Operations);
            U = obj.unsafeRegion;
            
            if obj.numCores > 1
                obj.start_pool;
                
            else
            
                S.data = obj.inputSet;
                S.opsIndex = 1; 
                safe = 'safe';
                while strcmp(safe, 'safe') && ~isempty(S)
                    S1 = S(1).data;
                    id = S(1).opsIndex;                    
                    if id < N
                        S2 = obj.Operations(id).execute(S1);
                        if length(S2) == 2
                            S3_1.data = S2(1);
                            S3_1.opsIndex = id + 1;
                            S3_2.data = S2(2);
                            S3_2.opsIndex = id + 1;
                            S3 = [S3_1 S3_2];
                            S(1) = [];
                            S = [S3 S];
                        else                            
                            S4.data = S2;
                            S4.opsIndex = id + 1;
                            S(1) = [];
                            S = [S4 S];
                        end
                    else
                        % checking safety of the leaf sets
                        S2 = obj.Operations(id).execute(S1);
                        if length(S2) == 2
                            H1 = S2(1).intersectHalfSpace(U.G, U.g);
                            H2 = S2(2).intersectHalfSpace(U.G, U.g);
                            if ~isempty(H1)
                                if strcmp(obj.reachMethod, 'exact-star')
                                    safe = 'unsafe';
                                    CEx = Star(obj.inputSet.V, H1.C, H1.d, H1.predicate_lb, H1.predicate_ub);
                                else
                                    safe = 'unknown';
                                    CEx = [];
                                end
                            elseif ~isempty(H2)
                                if strcmp(obj.reachMethod, 'exact-star')
                                    safe = 'unsafe';
                                    CEx = Star(obj.inputSet.V, H2.C, H2.d, H2.predicate_lb, H2.predicate_ub);
                                else
                                    safe = 'unknown';
                                    CEx = [];
                                end
                            end
                        else
                            H = S2.intersectHalfSpace(U.G, U.g);
                            if ~isempty(H)
                                if strcmp(obj.reachMethod, 'exact-star')
                                    safe = 'unsafe';
                                    CEx = Star(obj.inputSet.V, H.C, H.d, H.predicate_lb, H.predicate_ub);
                                else
                                    safe = 'unknown';
                                    CEx = [];
                                end
                            end
                        end
                        
                        S(1) = [];
                        
                    end
                    
                end
                
            end
            
            
            
        end
        
        
        
        
        
        
    end
    
    
    methods(Static) % parse Matlab trained feedforward network
        
        function nnvNet = parse(MatlabNet)
            % @MatlabNet: A feedforward network trained by Matlab
            % @nnvNet: an NNV FFNNS object for reachability analysis and
            % verification
            
            % author: Dung Tran
            % date: 9/13/2019
            
            if ~isa(MatlabNet, 'network')
                error('Input is not a matlab network');
            end
            
            n = length(MatlabNet.b); % number of layers
            
            if sum(MatlabNet.layerConnect, 'all') ~= n-1
                error('The network is not a feedfoward network');
            end
            
            Layers = [];
            for i=1:n
                act_func = MatlabNet.layers{i}.transferFcn;
                                              
                if i==1
                   
                   W = MatlabNet.IW{1};
                   b = MatlabNet.b{1};
                   Layer = LayerS(W, b, act_func);
                   
                else
                  
                   W = MatlabNet.LW{i, i-1};
                   b = MatlabNet.b{i};
                   Layer = LayerS(W, b, act_func);
                    
                end
                
                Layers = [Layers Layer];
                
            end
            
            nnvNet = FFNNS(Layers);
            
        end
        
    end
    
    
    
end

