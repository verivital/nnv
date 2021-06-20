classdef DAGNN < handle
    % DAG neural network class is a class for Verification of Directed
    % acyclic graph network
    % Author: Dung Tran
    % Date: 4/14/2020
    
    properties
        
        Name = 'dagnn'; % name of the network
        Layers = {}; % An array of Layers, eg, Layers = [L1 L2 ...Ln]
        Connections = []; % table of connection
        
        InputNames = [];
        OutputNames = [];
        
        numLayers = 0; % number of Layers
        InputSize = 0; % number of Inputs
        OutputSize = 0; % number of Outputs
        
        % properties for reachability analysis        
        reachMethod = 'approx-star';    % reachable set computation scheme, default - 'star'
        reachOption = []; % parallel option, default - non-parallel computing
        numCores = 0; % number of cores (workers) using in computation
        reachSet = [];  % reachable set for each layers
        reachTime = []; % computation time for each layers
        totalReachTime = 0; % total computation time
                
    end
    
    
    methods % constructor, evaluation, sampling, print methods
        
        % constructor
        function obj = DAGNN(varargin)
            
            switch nargin
                case 5
                    name = varargin{1};
                    layers = varargin{2};
                    connections = varargin{3};
                    inputNames = varargin{4};
                    outputNames = varargin{5};
                                        
                    % check input arguments
                    if ~ischar(name)
                        error('Invalid type of network name, should be char');
                    end
                    
                    if ~ismatrix(layers)
                        error('Layers should be an array of layers');
                    end
                    
                    if ~istable(connections)
                        error('Connection should be a table');
                    end
                    
                    if ~iscell(inputNames)
                        error('InputNames should be a cell');
                    end
                    
                    if ~iscell(outputNames)
                        error('OutputNames should be a cell');
                    end
                    
                    
                    nL = length(layers); % number of Layers                                       
                    obj.Name = name;
                    obj.Layers = layers;
                    obj.numLayers = nL;    % number of layers
                    obj.Connections = connections;
                    obj.InputNames = inputNames;
                    obj.OutputNames = outputNames;
                    if isprop(layers{1}, 'InputSize')
                        obj.InputSize = layers{1}.InputSize; % input size
                    end
                    if isprop(layers{nL}, 'OutputSize')
                        obj.OutputSize = layers{nL}.OutputSize; % output size
                    end
                
                case 3
                    name = varargin{1};
                    layers = varargin{2};
                    connections = varargin{3};
                    % check input arguments
                    if ~ischar(name)
                        error('Invalid type of network name, should be char');
                    end
                    
                    if ~ismatrix(layers)
                        error('Layers should be an array of layers');
                    end
                    
                    if ~istable(connections)
                        error('Connection should be a table');
                    end
                                       
                    
                    nL = length(layers); % number of Layers                                       
                    obj.Name = name;
                    obj.Layers = layers;
                    obj.numLayers = nL;    % number of layers
                    obj.Connections = connections;
                    if isprop(layers(1), 'InputSize')
                        obj.InputSize = layers(1).InputSize; % input size
                    end
                    if isprop(layers(nL), 'OutputSize')
                        obj.OutputSize = layers(nL).OutputSize; % output size
                    end
                
                case 2
                    
                    layers = varargin{1};
                    connections = varargin{2};
                    % check input arguments
                    if ~ismatrix(layers)
                        error('Layers should be an array of layers');
                    end
                    
                    if ~istable(connections)
                        error('Connection should be a table');
                    end
                                       
                    
                    nL = length(layers); % number of Layers                                       
                    obj.Layers = layers;
                    obj.numLayers = nL;    % number of layers
                    obj.Connections = connections;
                    if isprop(layers(1), 'InputSize')
                        obj.InputSize = layers(1).InputSize; % input size
                    end
                    if isprop(layers(nL), 'OutputSize')
                        obj.OutputSize = layers(nL).OutputSize; % output size
                    end
                
                    
                case 0
                    
                    obj.Name = 'emptyDAGN';
                    obj.Layers = [];
                    obj.numLayers = 0;    % number of layers
                    obj.Connections = [];
                    obj.InputNames = {};
                    obj.OutputNames = {};
                    
                otherwise
                    error('Invalid number of inputs, should be 0, 2, 3 or 5');
            end
                      
        end
                
        
        % Evaluation of a CNN
        function y = evaluate(varargin)
            % Evaluation of this FFNN
            % @x: input vector x
            % @y: output vector y
            % @layer_id: layer index
            
            switch nargin
                case 3
                    obj = varargin{1};
                    x = varargin{2};
                    layer_id = varargin{3};
                case 2
                    obj = varargin{1};
                    x = varargin{2};
                    layer_id = obj.numLayers;
                otherwise
                    error("Invalid number of input arguments, should be 1 or 2");
            end
            
            if layer_id < 1
                error('Invalid layer index');
            end
            
            y = x;
            for i=1:layer_id
                y = obj.Layers{i}.evaluate(y);
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
        
        function [IS, reachTime] = reach(varargin)            
            % @I: input set, a star set          
            % @method: = 'exact-star' or 'approx-star' -> compute reach set using stars
            %            'abs-dom' -> compute reach set using abstract
            %            domain (support in the future)
            %            'face-latice' -> compute reach set using
            %            face-latice (may support in the future)
            
            % @IS: output set is an ImageStar
            % @reachTime : reachable set computation time
            
            % author: Dung Tran
            % date:4/15/2020
            
            
            switch nargin 
                
                case 2
                    
                    obj = varargin{1};
                    inputSet = varargin{2};
                    obj.reachMethod = 'approx-star';
                    obj.numCores = 1; 
                    
                case 3 
                    
                    obj = varargin{1};
                    inputSet = varargin{2};
                    obj.reachMethod = varargin{3};
                    obj.numCores = 1; 
                    
                case 4
                    
                    obj = varargin{1};
                    inputSet = varargin{2};
                    obj.reachMethod = varargin{3};
                    obj.numCores = varargin{4}; 
                 
                otherwise 
                    
                    error('Invalid number of input arguments, the number should be 1, 2 or 3');
                
            end       
            
            if  obj.numCores > 1
                obj.start_pool;
                obj.reachOption = 'parallel';
            else
                obj.reachOption = [];
            end
            
            obj.reachSet = cell(1, obj.numLayers+1);
            obj.reachTime = zeros(1, obj.numLayers);
            fprintf('\nPerform reachability analysis for the network %s...', obj.Name);
            obj.reachSet{1} = inputSet;
            for i=2:obj.numLayers+1
                fprintf('\nPerforming analysis for Layer %d (%s)...', i-1, obj.Layers{i-1}.Name);
                start_time = tic;
                obj.reachSet{i} = obj.Layers{i-1}.reach(obj.reachSet{i-1}, obj.reachMethod, obj.reachOption);
                obj.reachTime(i-1) = toc(start_time);
                fprintf('\nReachability analysis for Layer %d (%s) is done in %.5f seconds', i-1, obj.Layers{i-1}.Name, obj.reachTime(i-1));
                fprintf('\nThe number of reachable sets at Layer %d (%s) is: %d', i-1, obj.Layers{i-1}.Name, length(obj.reachSet{i}));
            end
            fprintf('\nReachability analysis for the network %s is done in %.5f seconds', obj.Name, sum(obj.reachTime));
            fprintf('\nThe number ImageStar in the output sets is: %d', length(obj.reachSet{obj.numLayers+1}));
            obj.totalReachTime = sum(obj.reachTime);
            IS = obj.reachSet{obj.numLayers+1};
            reachTime = obj.totalReachTime;
        end
        
        
    end
    
    
    
    
    methods(Static)
       
        % parse a dag neural network from matlab for reachability analysis
        function dagnn = parse(varargin)
            % @net: input network, should be a dagnetwork or lgraph object
            % @dagnn: the dagnn network for reachability analysis
            
            % the constructed DAGNN for reachability analysis get rid of the
            % these following layers:
            % 2) Dropout Layer (is not used for prediction phase)
            % 3) Softmax Layer
            
            % author: Dung Tran
            % date: 4/14/2020
            % update: 4/10/2020
            
            switch nargin
                case 1
                    net = varargin{1};
                    name = 'parsed_net';
                case 2
                    net = varargin{1};
                    name = varargin{2};
                otherwise
                    error('Invalid number of input arguments, should be 1 or 2');
            end
            
            
            
            if ~isa(net, 'DAGNetwork') && ~isa(net, 'LayerGraph')
                error('the parsed object is not a Matlab DAGNetwork object or a LayerGraph object');
            end
            
            n = length(net.Layers); % number of layers
                        
            Ls = {};
            j = 0;
            for i=1:n
                L = net.Layers(i);
                fprintf('\nParsing Layer %d...', i);
                
                if isa(L, 'nnet.cnn.layer.DropoutLayer') || isa(L, 'nnet.cnn.layer.ClassificationOutputLayer')                   
                    fprintf('\nLayer %d is a %s class which is neglected in the analysis phase', i, class(L));
                else
                    
                    if isa(L, 'nnet.cnn.layer.ImageInputLayer')
                        Li = ImageInputLayer.parse(L);
                    elseif isa(L, 'nnet.cnn.layer.Convolution2DLayer') 
                        Li = Conv2DLayer.parse(L);
                    elseif isa(L, 'nnet.cnn.layer.ReLULayer')
                        Li = ReluLayer.parse(L);
                    elseif isa(L, 'nnet.cnn.layer.BatchNormalizationLayer')
                        Li = BatchNormalizationLayer.parse(L);
                    elseif isa(L, 'nnet.cnn.layer.MaxPooling2DLayer')
                        Li = MaxPooling2DLayer.parse(L);
                    elseif isa(L, 'nnet.cnn.layer.AveragePooling2DLayer')
                        Li = AveragePooling2DLayer.parse(L);
                    elseif isa(L, 'nnet.cnn.layer.FullyConnectedLayer')
                        Li = FullyConnectedLayer.parse(L);
                    elseif isa(L, 'nnet.cnn.layer.MaxUnpooling2DLayer')
                        Li = MaxUnpooling2DLayer.parse(L);
                    elseif isa(L, 'nnet.cnn.layer.PixelClassificationLayer')
                        Li = PixelClassificationLayer.parse(L);
                    elseif isa(L, 'nnet.cnn.layer.SoftmaxLayer')
                        Li = SoftmaxLayer.parse(L);
                    else                     
                        fprintf('\nLayer %d is a %s which have not supported yet in nnv, please consider removing this layer for the analysis', i, class(L));
                        error('\nUnsupported Class of Layer');                     
                    end
                    j = j + 1;
                    Ls{j} = Li;       
                end
                            
            end
            
            % update connection table                           
            dagnn = DAGNN(name, Ls, net.Connections, net.InputNames, net.OutputNames);
            fprintf('\nParsing network is done successfully and %d Layers are neglected in the analysis phase', n - j);
            
        end
        
             
        %
        
        
    end
    
    
    
end

