 classdef FFNN < handle
    % FeedForward Neural Network Class
    % Dung Tran: 8/22/2018
    
    properties
        Layers = []; % An array of Layers, eg, Layers = [L1 L2 ...Ln]
        nL = 0; % number of Layers
        nN = 0; % number of Neurons
        nI = 0; % number of Inputs
        nO = 0; % number of Outputs
        
        % properties for reach set computation
        
        reachScheme = '';    % reachable set computation scheme
        numOfCores = 0; % number of cores (workers) using in computation
        inputSet = [];  % input set
        reachSet = [];  % reachable set for each layers
        outputSet = []; % output reach set
        reachTime = []; % computation time for each layers
        totalReachTime = 0; % total computation time
        n_ReLU_reduced = []; % number of ReLU stepReach operation is reduced
        total_n_ReLU_reduced = 0; % total number of ReLU stepReach operation is reduced
        
    end
    
    methods
        
        % constructor
        function obj = FFNN(Layers)
            nL = size(Layers, 2); % number of Layer
            for i=1:nL
                L = Layers(i);
                if ~isa(L, 'Layer')
                    error('Element %d of Layers array is not a Layer object', i);
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
        
        % Reach Set computation
        function [R, t] = reach(obj, I, scheme, numOfCores, max_nP)
            % Reach Set Computation of this FFNN
            
            % @I: input set which is a polyhedron
            
            % @scheme: = 'exact' -> compute the exact reach set
            %          - if I are Stars, we use exact method with star set
            %          - if I are Polyhedron, we use exact method with
            %          polyhedron set.
            %          = 'approx' -> compute the over-approximate reach
            %          set using boxes
            %          = 'mix' -> compute the over-approximate reach set
            %          using mixing scheme, exact-approx
            
            % @R: output set which is an array of polyhedron if option = 'exact'
            %     or is a polyhedron if option = 'approx'
            
            % @rn: number of ReLU_i (stepReach) operation reduced
            
            % @t : computation time which is an array containing computation time
            %     on each layer.
            
            % @numOfCores: number of cores you want to run the reachable
            % set computation, @numOfCores >= 1, maximum is the number of
            % cores in your computer. We also can run this method on a
            % clusters. You need to change: parpool('local', numOfCores) to
            % parcluster(your_cluster_profile), you also need to set up
            % your local clusters with installed MPT toolbox also. see: 
            % https://www.mathworks.com/help/distcomp/discover-clusters-and-use-cluster-profiles.html
           
            % @max_nP: the maximum number of Polyhedra or boxes used to
            % over-approximate the reachable set of each layers
            % this input is only used for 'approx' and 'mix' schemes.
           
            
            
        
            if numOfCores < 1
                error('Number of cores should be at least one');
            end
            
            if numOfCores == 1
                parallel = 'single';
            else
                parallel = 'parallel';
            end
            
            obj.numOfCores = numOfCores;
            obj.reachScheme = scheme;
            obj.inputSet = I;
            obj.reachSet = cell(1, obj.nL);
            obj.reachTime = []; 
            
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
            end            
            
            In = I;
                        
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
                    if isnan(estimatedTime)
                        estimatedTime = 0.0;
                    end
                    fprintf('\nEstimated computation time: ~ %.5f seconds', estimatedTime);
                end
                
                if strcmp(scheme, 'exact')                   

                    [In, rn1, t1] = obj.Layers(i).reach_exact(In, parallel);
                    obj.n_ReLU_reduced = [obj.n_ReLU_reduced rn1];                        
                    
                else
                                        
                    if strcmp(scheme, 'approx')

                        [In, t1] = obj.Layers(i).reach_approx_box_coarse(In, max_nP, parallel);

                    elseif strcmp(scheme, 'mix')

                        [In, t1] = obj.Layers(i).reach_mix(In, max_nP, parallel);
                        
                    %elseif strcmp(scheme, 'approx-star')

                    %    [In, rn1, t1] = obj.Layers(i).reach_approx_star(In, parallel);  
                    %    obj.n_ReLU_reduced = [obj.n_ReLU_reduced rn1];
                        
                    else      
                        error('Unknown scheme');
                    end
                end

                obj.reachSet{1, i} = In;
                obj.reachTime = [obj.reachTime t1];
                
                fprintf('\nExact computation time: %.5f seconds', t1);               
                fprintf('\nNumber of polyhedra at the output of layer %d: %d', i, length(In));
                
                               
            end
            
            obj.outputSet = In;
            obj.totalReachTime = sum(obj.reachTime);
            obj.total_n_ReLU_reduced = sum(obj.n_ReLU_reduced);
            R = obj.outputSet;
            t = obj.totalReachTime;
                       
            fprintf('\nTotal reach set computation time: %.5f seconds', obj.totalReachTime);
            if ~isempty(obj.n_ReLU_reduced)
                fprintf('\nTotal number of ReLU_i stepReach operations reduced: %d', obj.total_n_ReLU_reduced);
            end
            fprintf('\nTotal number of polyhedron at the output layer: %d', length(obj.outputSet));
        end
        
        
        % reachability analysis with conservativeness guarantee,
        % this is used for lazy-approximate + input partition scheme
        % and mixing scheme
        % NOTE : *** This method only works for hyper-rectangle input set 
        function [R, t] = reach_approx_with_CSV_guarantee(obj, I, desired_csv, k_max, numOfCores, n_samples)
            % @I: input set
            % @desired_csv: desired conservativeness specified by user
            % the algorithm stop when all outputs ranges's conservativeness <= desired_csv
            % @k_max: maximum divide times, for example k=1 -> each x[i]
            % with be divide into 2 subsegments, k=2 -> each x[i] is
            % divided into 4 subsegments, and so on.
            % @numOfCores: number of cores used in computation
            % @n_samples: number of samples used to estimate the unknown
            % actual ranges of the outputs
            
            % @R: output reachable set with desired conservativeness
            % @t: computation time
            
            t1 = tic;
            obj.reach(I, 'approx', numOfCores, []);
            csv_vec = obj.estimate_CSV(I, n_samples);
            k = 0;
            while (max(csv_vec) > desired_csv) && k < k_max
                k = k + 1;
                fprintf('\nPartitioning input space...')
                I1 = Partition.partition_box(I, k);
                fprintf('\nRecompute reachable set with %d smaller input sets', length(I1));
                obj.reach(I1, 'approx', numOfCores, []);
                csv_vec = obj.estimate_CSV(I, n_samples);
            end
            if k >= k_max
                fprintf('\nParitioning reaches to the limitation');
                if max(csv_vec) > desired_csv
                    fprintf('\nThe Reachable set refinement is not successful, try again by increasing the dividing limitation k_max');
                    fprintf('\nThe conservativeness (in percentage) of the return output reachable set is:');
                    display(csv_vec);
                else
                    fprintf('\nReachable set refinement is done successfully');
                    fprintf('\nThe conservativeness (in percentage) of the return output reachable set is:');
                    display(csv_vec);
                end
            else
                fprintf('\nReachable set refinement is done successfully');
                fprintf('\nThe conservativeness (in percentage) of the return output reachable set is:');
                display(csv_vec);
            end
            R = obj.outputSet;
            t = toc(t1);
            
        end
        
        % estimate conservativeness of output reachable set compared with
        % the actual output ranges
        function [csv_vec, r, computed_range, est_range] = estimate_CSV(obj, I, n_samples)
            % @I: input set, currently only support input set as a box
            % @n_samples: number of samples used to estimate the actual
            % ranges 
            % @CSV_vec: conservativeness vector
            % @r: = 0: exact range
            %     = 1: over-approximation
            %     = 2: under-approximation
            %     = 3: wrong range
            % @est_range: estimated range from sampling the network
            % @computed_range: the computed range from reachable set
            % computation
                        
            if isempty(obj.outputSet)
                error('output set is empty, call reach method to compute the reachable set first');
            end
               
            B = Reduction.hypercubeHull(obj.outputSet);
            computed_range = [B.lb B.ub];
            est_range = obj.estimate_ranges(I, n_samples);
            [csv_vec, r] = CSV.getConservativeness(computed_range, est_range);     
            
        end
        
        
        % Sample of FFNN
        function Y = sample(obj, V)
            % sample the output of each layer in the FFNN based on the
            % vertices of input set I, this is useful for testing.
            % @V : array of vertices to evaluate
            % @Y : output which is a cell array
        
            Y = cell(1, obj.nL);
            In = V;
            for i=1:obj.nL
                In = obj.Layers(i).sample(In);
                Y{1, i} = In;
            end
        end
        
        % quick estimate ranges of an FFNN
        function est_ranges = estimate_ranges(obj, I, n_samples)
            % @I: input set, currently works for a hyper-rectangle
            % @n_samples: number of samples used to estimate the output
            % ranges
            
            % @est_range: the estimated ranges of the outputs with given
            % input set I.
            
            if ~isa(I, 'Box')
                I.outerApprox;
                lb = I.Internal.lb;
                ub = I.Internal.ub;
            else
                lb = I.lb;
                ub = I.ub;
            end
            
            V1 = Reduction.getVertices(lb, ub); % get vertices from the set,
            % these vertices usually are the input points corresponding to the 
            % maximum and minimum values of the output ranges.
            
            % sampling the network with n_samples of input vector           
            % get sampled input vectors
            X = cell(1, obj.nI);
            V = [];
            for i=1:obj.nI
                X{1, i} = (ub(i) - lb(i)).*rand(n_samples, 1) + lb(i);
                V = vertcat(V, X{1, i}');
            end
      
            V = horzcat(V, V1);
                        
            Y = obj.sample(V); % evaluate the network with a set of input vectors V
            output = Y{1, obj.nL}; % sampled output vectors
            output_lb = min(output, [], 2); % estimate ranges of outputs
            output_ub = max(output, [], 2);
            
            est_ranges = [output_lb output_ub];
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
                fprintf(f, '\nReachability scheme: %s', obj.reachScheme);
                fprintf(f, '\nNumber of cores used in computation: %d', obj.numOfCores);
                
                for i=1:length(obj.reachSet)-1
                    fprintf(f, '\nLayer %d reach set consists of %d polyhedrons that are computed in %.5f seconds', i, length(obj.reachSet{1, i}), obj.reachTime(i));
                    if ~isempty(obj.n_ReLU_reduced)
                        fprintf(f, '\nNumber of ReLU_i stepReach operations neglected for layer %d: %d', i, obj.n_ReLU_reduced(i));
                    end
                end
                fprintf(f, '\nOutput Layer reach set consists of %d polyhedrons that are computed in %.5f seconds', length(obj.reachSet{1, obj.nL}), obj.reachTime(obj.nL)); 
                fprintf(f, '\nTotal reachable set computation time: %.5f', obj.totalReachTime);
                
                if ~isempty(obj.n_ReLU_reduced)
                    fprintf(f, '\nTotal number of ReLU_i stepReach operations neglected: %d', obj.total_n_ReLU_reduced);
                end
                
            end
               
        end
        
        
    end
    
end

