classdef FFNNS < handle
    % FFNNS Class is a new feedforward network class used to replace the old
    % FFNN in the future, FFNNS does not support polyhedron-based
    % reachability analysis methods.
    % Author: Dung Tran
    % Date: 27/2/2019
    
    properties
        Layers = []; % An array of Layers, eg, Layers = [L1 L2 ...Ln]
        nL = 0; % number of Layers
        nN = 0; % number of Neurons
        nI = 0; % number of Inputs
        nO = 0; % number of Outputs
        
        % properties for reach set computation
        
        reachMethod = 'star';    % reachable set computation scheme, default - 'star'
        reachOption = []; % parallel option, default - non-parallel computing
        numCores = 0; % number of cores (workers) using in computation
        inputSet = [];  % input set
        reachSet = [];  % reachable set for each layers
        outputSet = []; % output reach set
        reachTime = []; % computation time for each layers
        totalReachTime = 0; % total computation time
        
    end
    
    
    methods % constructor, evaluation, sampling, print methods
        
        % constructor
        function obj = FFNNS(Layers)
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
                fprintf(f, '\nReachability method: %s', obj.reachMethod);
                fprintf(f, '\nNumber of cores used in computation: %d', obj.numCores);
                
                for i=1:length(obj.reachSet)-1
                    fprintf(f, '\nLayer %d reach set consists of %d sets that are computed in %.5f seconds', i, length(obj.reachSet{1, i}), obj.reachTime(i));
                end
                fprintf(f, '\nOutput Layer reach set consists of %d sets that are computed in %.5f seconds', length(obj.reachSet{1, obj.nL}), obj.reachTime(obj.nL)); 
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
            % @I: input set          
            % @method: = 'star' -> compute reach set using stars
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
                    obj.numCores = varargin{3};
                case 2
                    obj = varargin{1};
                    obj.inputSet = varargin{2};
                    obj.numCores = 1;
                    
                otherwise
                    error('Invalid number of input arguments (should be 2 or 4)');
            end
            
            
            if obj.numCores < 1
                error('Number of cores should be at least one');
            elseif obj.numCores == 1
                obj.reachOption = []; % don't use parallel computing
            else
                obj.reachOption = 'parallel';
            end
            
            obj.reachSet = cell(1, obj.nL);
            obj.reachTime = [];
            
            % start parallel pool in Matlab
            if obj.numCores > 1
                obj.start_pool;
            end
                        
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
                
                obj.reachSet{1, i} = In;
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
    
    
end

