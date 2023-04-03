classdef BNN < handle
    % BNN Class is a new binary neural network
    % reachability analysis method: not implemented
    %
    % author: Michael Ivashchenko
    
    properties
        Name = 'net';
        Layers = {}; % An array of Layers, eg, Layers = {L1 L2 ... Ln}
        nL = 0; % number of Layers
        nN = 0; % number of Neurons
        nP = 0; % number of Parameters
        nI = 0; % number of Inputs
        nO = 0; % number of Outputs
        
        % properties for reach set computation
        
        reachMethod = 'exact-star';    % reachable set computation scheme, default - 'star'
        reachOption = []; % parallel option, default - non-parallel computing
        relaxFactor = 0; % use only for approximate star method, 0 mean no relaxation
        numCores = 1; % number of cores (workers) using in computation
        inputSet = [];  % input set
        reachSet = [];  % reachable set for each layer
        outputSet = []; % output reach set
        reachTime = []; % computation time for each layer
        signReachTime = []; % computation time for each sign layer
        numReachSet = []; % number of reach sets over layers
        numSamples = 0; % default number of samples using to falsify a property
        unsafeRegion = []; % unsafe region of the network
        getCounterExs = 0; % default, not getting counterexamples
        
        signLayerReachTime = 0;
        signLayerCounter = 0;
        
        Operations = []; % flatten a network into a sequence of operations
        
        dis_opt = []; % display option
        lp_solver = 'linprog'; % lp solver option, should be glpk or linprog
        
        type = 'mlp';
        
        supportedLayers = [
               "LayerS"; "BatchNormalizationLayer"; ...
               "Conv2DLayer"; "MaxPooling2DLayer"; "FlattenLayer"; ...
               "SignLayer"; "FullyConnectedLayer"
            ];
    end
    
    methods % constructor, evaluation
        
        % constructor
        function obj = BNN(varargin)
            switch nargin
                case 0
                    obj.Name = 'bnn';
                    return
                case 1
                    Layers = varargin{1};
                    obj.Name = 'bnn';
                case 2
                    Layers = varargin{1};
                    obj.Name = 'bnn';
                    obj.dis_opt = [obj.dis_opt; varargin{2}];
                otherwise
                    error('Invalid number of inputs');
            end
            
            % set up the number of Layers
            nL = length(Layers);
            if any(strcmp(obj.dis_opt(:), 'display'))
                disp(["Layers number: " num2str(nL)])
            end
            
            % validate Layers types
            isSupported = false;
            for i=1:nL               
                for j=1:length(obj.supportedLayers)
                    if isa(Layers{i}, obj.supportedLayers(j))
                        isSupported = true;
                        break;
                    end
                end
                
                if ~isSupported
                    error('Element %d of Layers array is not supported', i)
                end
                
                isSupported = false;
            end
            if any(strcmp(obj.dis_opt(:), 'display'))
                disp("Layers type check successful")
            end
            
            obj.Layers = Layers; % set of Layers
            obj.nL = nL; % number of Layers
            
            for i=1:nL
                if ~isa(Layers{i}, "BatchNormalizationLayer") && ~isa(Layers{i}, "FlattenLayer") && ~isa(Layers{i}, "SignLayer") && obj.nI == 0                    
                    
                    if isa(Layers{i}, "LayerS") || isa(Layers{i}, "FullyConnectedLayer")
                        parSize = size(Layers{i}.W);
                        obj.nI = size(Layers{i}.W, 2);
                    elseif isa(Layers{i}, "Conv2DLayer")
                        parSize = size(Layers{i}.Weights);
                        obj.nI = size(Layers{i}.Weights(:,1:3));
                    end
                    
                    
                    obj.nP = obj.nP + prod(parSize, 'all');
                    if sum(Layers{i}.b) ~= 0
                        bSize = size(Layers{i}.b);
                        obj.nP = obj.nP + prod(bSize, 'all');
                    end
                end
                
                if ~isa(Layers{nL - i + 1}, "BatchNormalizationLayer") && ~isa(Layers{nL - i + 1}, "FlattenLayer") && ~isa(Layers{nL - i + 1}, "SignLayer") && obj.nO == 0
                    parSize = size(Layers{nL - i + 1}.b);
                    obj.nO = prod(parSize, 'all');
                end
            end
        end
        
        % Evaluation of a BNN
        function y = evaluate(obj, x)
            % Evaluates current BNN
            % @x: input vector x
            % @y: output vector y
            
            y = x;
            for i=1:obj.nL
                if any(strcmp(obj.dis_opt(:), 'display'))
                    disp(["Evaluating layer " i]);
                    disp(obj.Layers{i});
                end
                
                if isa(obj.Layers{i}, 'Conv2DLayer')
                    y = obj.Layers{i}.evaluate2(y, 'double');
                else
                    y = obj.Layers{i}.evaluate(y);
                    if (isa(obj.Layers{i}, 'BatchNormalizationLayer') && length(size(y)) == 3 && size(y, 1) == 1 && size(y,2) == 1) %|| isa(obj.Layers{i}, 'FlattenLayer'))  
                        y = reshape(y, [size(y, 3) 1]);
                    end
                end
            end
        end
        
        % reachability method
        function [S, bnnReachTime, signReachTime] = reach(varargin)            
            % @I: input set, a star set          
            % @method: = 'exact-star' or 'approx-star' -> compute reach set using stars
            %            'approx-zono' -> compute reach set using zonotope
            %            domain (support in the future)
            
            % @R: output set 
            % @t : computation time 
            
            % @numOfCores: number of cores you want to run the reachable
            % set computation, @numOfCores >= 1, maximum is the number of
            % cores in your computer. We also can run this method on a
            % clusters. You need to change: parpool('local', numOfCores) to
            % parcluster(your_cluster_profile), you also need to set up
            % your local clusters with installed MPT toolbox also. see: 
            % https://www.mathworks.com/help/distcomp/discover-clusters-and-use-cluster-profiles.html
           
            % author: Mykhailo Ivashchenko
            
            % parse inputs 
            switch nargin
                
                 case 7
                    obj = varargin{1}; % FFNNS object
                    obj.inputSet = varargin{2}; % input set
                    obj.reachMethod = varargin{3}; % reachability analysis method
                    obj.numCores = varargin{4}; % number of cores used in computation
                    obj.relaxFactor = varargin{5}; 
                    obj.dis_opt = varargin{6};
                    obj.lp_solver = varargin{7};
          
                case 6
                    obj = varargin{1}; % FFNNS object
                    obj.inputSet = varargin{2}; % input set
                    obj.reachMethod = varargin{3}; % reachability analysis method
                    obj.numCores = varargin{4}; % number of cores used in computation
                    obj.relaxFactor = varargin{5}; 
                    obj.dis_opt = varargin{6}; 
                
                case 5
                    obj = varargin{1}; % FFNNS object
                    obj.inputSet = varargin{2}; % input set
                    obj.reachMethod = varargin{3}; % reachability analysis method
                    obj.numCores = varargin{4}; % number of cores used in computation
                    obj.relaxFactor = varargin{5}; % used only for approx-star method
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
                    if ~isstruct(varargin{2})
                        obj.inputSet = varargin{2};
                        obj.reachMethod = 'exact-star';
                        obj.numCores = 1;
                    else
                        if isfield(varargin{2}, 'inputSet')
                            obj.inputSet = varargin{2}.inputSet;
                        end
                        if isfield(varargin{2}, 'numCores')
                            obj.numCores = varargin{2}.numCores;
                        end
                        if isfield(varargin{2}, 'reachMethod')
                            obj.reachMethod = varargin{2}.reachMethod;
                        end
                        if isfield(varargin{2}, 'dis_opt')
                            obj.dis_opt = varargin{2}.dis_opt;
                        end
                        if isfield(varargin{2}, 'relaxFactor')
                            obj.relaxFactor = varargin{2}.relaxFactor;
                        end
                        if isfield(varargin{2}, 'lp_solver')
                            obj.lp_solver = varargin{2}.lp_solver;
                        end
                        
                    end   
                        
                    
                otherwise
                    error('Invalid number of input arguments (should be 1, 2, 3, 4, 5, or 6)');
            end
            
            obj.outputSet = [];
            obj.reachTime = [];
            obj.signReachTime = [];
            
            % if reachability analysis method is an over-approximate
            % method, we use 1 core for computation
            if strcmp(obj.reachMethod, 'approx-star')
                obj.numCores = 1;
            end
            
            if obj.numCores == 1
                obj.reachOption = []; % don't use parallel computing
            else
                obj.start_pool;       % start parallel pool in Matlab
                obj.reachOption = 'parallel';
            end
            
            obj.reachSet = cell(1, obj.nL);
            obj.numReachSet = zeros(1, obj.nL);
                        
            % compute reachable set
            In = obj.inputSet;
                        
            for i=1:obj.nL
                if any(strcmp(obj.dis_opt(:), 'display'))
                    fprintf('\nComputing reach set for Layer %d ...\n', i);
                    disp(obj.Layers{i});
                end
                
                st = tic;
                
                if isa(obj.Layers{i}, 'LayerS') && strcmp(obj.type, 'xnor') && isa(In, 'ImageStar')
                   In = In.toStar();
                end
                
                
                if (isa(obj.Layers{i}, 'SignLayer') && strcmp(obj.type, 'xnor'))
                    In = obj.Layers{i}.reach(In, 'approx-star', obj.reachOption, obj.relaxFactor, obj.dis_opt, obj.lp_solver);
                else
                    In = obj.Layers{i}.reach(In, obj.reachMethod, obj.reachOption, obj.relaxFactor, obj.dis_opt, obj.lp_solver);
                end
                
                t1 = toc(st);
                
                if isa(obj.Layers{i}, 'SignLayer')
                    obj.signReachTime = [obj.signReachTime t1];
                end
                
                obj.numReachSet(i) = length(In);
                obj.reachTime = [obj.reachTime t1];
                
                if any(strcmp(obj.dis_opt(:), 'display'))
                    fprintf('\nReachable set computation time: %.5f seconds\n', t1);               
                    fprintf('\nNumber of reachable set at the output of layer %d: %d\n', i, length(In));               
                end           
            end
            
            obj.outputSet = In;
            
            S = obj.outputSet;
            bnnReachTime = sum(obj.reachTime);
            signReachTime = sum(obj.signReachTime);
            
            if any(strcmp(obj.dis_opt(:), 'display'))
                fprintf('\nTotal reach set computation time: %.5f seconds\n', sum(obj.reachTime));
                fprintf('\nTotal sign reach set computation time: %.5f seconds\n', sum(obj.signReachTime));
                fprintf('\nAvg sign reach set computation time: %.5f seconds\n', mean(obj.signReachTime));
                fprintf('\nTotal number of output reach sets: %d\n', length(obj.outputSet));
            end
            
        end
        
        % start parallel pool for computing
        function start_pool(varargin)
            
            switch nargin
                case 1
                    obj = varargin{1};
                    nCores = obj.numCores;
                case 2
                    nCores = varargin{2};
                otherwise
                    error('Invalid number of input arguments');
            
            end
            
            if nCores > 1
                poolobj = gcp('nocreate'); % If no pool, do not create new one.
                if isempty(poolobj)
                    parpool('local', nCores); 
                else
                    if poolobj.NumWorkers ~= nCores
                        delete(poolobj); % delete the old poolobj
                        parpool('local', nCores); % start the new one with new number of cores
                    end                    
                end
            end   
            
        end
                
        % Prints the contents of the network
        function print(obj, file_name)
            % @file_name: name of file the network info is stored in
            
            file = fopen(file_name, 'w');
            fprintf(file, 'Binary Neural Network Information\n');
            fprintf(f, '\nNumber of layers: %d', obj.nL);
            fprintf(f, '\nNumber of neurons: %d', obj.nN);
            fprintf(f, '\nNumber of inputs: %d', obj.nI);
            fprintf(f, '\nNumber of outputs: %d', obj.nO);
        end
    end
end