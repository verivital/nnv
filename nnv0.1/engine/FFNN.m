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
        function [R, t] = reach(obj, I, scheme, numOfCores, numOfPolyhedra)
            % Reach Set Computation of this FFNN
            
            % @I: input set which is a polyhedron
            
            % @scheme: = 'exact' -> compute the exact reach set
            %          = 'approx-box' -> compute the over-approximate reach
            %          set using boxes
            %          = 'approx-polyhedron' -> compute the
            %          over-approximate reach set using polyhedra.
            
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
           
            % @numOfPolyhedra: the number of Polyhedra or boxes used to
            % over-approximate the reachable set of each layers
            % this input is only used for over-approximation scheme i.e.,
            % 'approx-box' or 'approx-polyhedron'
            %  @numOfPolyhedra is an array, for example [3 2 1 5] means, 
            % we use 3 boxes for over-approximating the reachable set of
            % layer 1, and 2 boxes for layer 2 ....
            
            
        
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
            
            % set up parallel computing with number of cores (workers)
            if numOfCores > 1
                poolobj = gcp('nocreate'); % If no pool, do not create new one.
                if isempty(poolobj)
                    poolobj = parpool('local', numOfCores); 
                else
                    if poolobj.NumWorkers ~= numOfCores
                        delete(poolobj); % delete the old poolobj
                        poolobj = parpool('local', numOfCores); % start the new one with new number of cores
                    end                    
                end
            end
            
            
            % check consistency of numOfPolyhedra if the scheme is an
            % over-approximate scheme.
            if ~strcmp(scheme, 'exact')

                [n1, m1] = size(numOfPolyhedra);
                if n1 ~=1
                    error('NumOfPolyhedra shoule be a one row array');
                end

                if m1 ~= obj.nL
                    error('Number of column of NumOfPolyhedra array should be equal the number of layers')
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
                    
                    fprintf('\nEstimated computation time: ~ %.5f seconds', estimatedTime);
                end
                
                if strcmp(scheme, 'exact')
                    
                    [In, rn1, t1] = obj.Layers(i).reach_exact(In, parallel);
                    obj.n_ReLU_reduced = [obj.n_ReLU_reduced rn1];
                    
                else
                                        
                    if numOfPolyhedra(i) < 1
                        error('Number of Polyhedra using to over-approximating reach set should be >= 1');
                    end
                    
                    if strcmp(scheme, 'approx-box')

                        %[In, t1] = obj.Layers(i).reach_approx_box(In, 1000, parallel);
                        [In, t1] = obj.Layers(i).reach_mix(In, 800, parallel);

                    elseif strcmp(scheme, 'approx-polyhedron')

                        [In, t1] = obj.Layers(i).reach_approx_polyhedron(In, numOfPolyhedra(i), 10, parallel);

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

