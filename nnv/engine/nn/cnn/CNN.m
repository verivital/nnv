classdef CNN < handle
    % CNN Class is a class for Verification of Convolutional Neural
    % Networks
    % Author: Dung Tran
    % Date: 6/27/2019
    
    properties
        
        Name = 'cnn'; % name of the network
        Layers = {}; % An array of Layers, eg, Layers = [L1 L2 ...Ln]
        numLayers = 0; % number of Layers
        numNeurons = 0; % number of Neurons
        InputSize = 0; % number of Inputs
        OutputSize = 0; % number of Outputs
        
        % properties for reach set computation
        
        reachMethod = 'approx-star';    % reachable set computation scheme, default - 'star'
        reachOption = []; % parallel option, default - non-parallel computing
        numCores = 0; % number of cores (workers) using in computation
        reachSet = [];  % reachable set for each layers
        reachTime = []; % computation time for each layers
        totalReachTime = 0; % total computation time
        
        features = {}; % outputs of each layer in an evaluation
        
    end
    
    
    methods % constructor, evaluation, sampling, print methods
        
        % constructor
        function obj = CNN(varargin)
            
            switch nargin
                case 4
                    name = varargin{1};
                    Layers = varargin{2};
                    inputsize = varargin{3};
                    outputsize = varargin{4};
                    nL = length(Layers); % number of Layers
                    for i=1:nL
                        L = Layers{i};
                        if ~isa(L, 'ImageInputLayer') && ~isa(L, 'AveragePooling2DLayer') && ~isa(L, 'Conv2DLayer') && ~isa(L, 'FullyConnectedLayer') && ~isa(L, 'MaxPooling2DLayer') && ~isa(L, 'ReluLayer')
                            fprintf('\nCurrent version of NNV supports ImageInputLayer, AveragePooling2DLayer, Convolutional2DLayer, FullyConnectedLayer, MaxPooling2DLayer, AveragePooling2DLayer and ReluLayer');
                            error('Element %d of Layers is not among supported layers in NNV', i);
                        end
                    end
                    
                    obj.Name = name;
                    obj.Layers = Layers;
                    obj.numLayers = nL;    % number of layers
                    obj.InputSize = inputsize; % input size
                    obj.OutputSize = outputsize; % output size
                    
                case 0
                    
                    obj.Layers = {};
                    obj.numLayers = 0;
                    obj.InputSize = 0;
                    obj.OutputSize = 0;
                    
                otherwise
                    error('Invalid number of inputs, should be 0 or 4');
            end
                      
        end
                
        
        % Evaluation of a CNN
        function [y, features] = evaluate(obj, x)
            % Evaluation of this FFNN
            % @x: input vector x
            % @y: output vector y
            % @features: output of all layers
            
            y = x;
            for i=1:obj.numLayers
                y = obj.Layers{i}.evaluate(y);
                obj.features{i} = y;
            end
            features = obj.features;      
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
            %            face-latice (support in the future)
            
            % @IS: output set is an ImageStar
            % @reachTime : reachable set computation time
            
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
       
        % parse a network from matlab for reachability analysis
        function cnn = parse(net, name)
            % @net: input network
            % @cnn: the cnn network for reachability analysis
            
            % the constructed CNN for reachability analysis get rid of the
            % these following layers:
            % 1) InputImageLayer
            % 2) Dropout Layer (is not used for prediction phase)
            % 3) Softmax Layer
            % 4) Classification Output Layer         
            
            % author: Dung Tran
            % date: 7/16/2019
            
            
            n = length(net.Layers); % number of layers
                                   
            Ls = {};
            inputSize = [];
            outputSize = [];
            
            j = 0; % counter of number of layers
            for i=1:n
                L = net.Layers(i);
                if isa(L, 'nnet.cnn.layer.DropoutLayer') || isa(L, 'nnet.cnn.layer.SoftmaxLayer') || isa(L, 'nnet.cnn.layer.ClassificationOutputLayer')                  
                    fprintf('\nLayer %d is a %s class which is neglected in the analysis phase', i, class(L));                   
                    if isa(L, 'nnet.cnn.layer.ImageInputLayer')
                        inputSize = L.InputSize;
                    elseif isa(L, 'nnet.cnn.layer.ClassificationOutputLayer')
                        outputSize = L.OutputSize;
                    end
                    
                else
                    
                    fprintf('\nParsing Layer %d...', i);
                    
                    if isa(L, 'nnet.cnn.layer.ImageInputLayer')
                        Li = ImageInputLayer.parse(L);
                    elseif isa(L, 'nnet.cnn.layer.Convolution2DLayer') 
                        Li = Conv2DLayer.parse(L);
                    elseif isa(L, 'nnet.cnn.layer.ReLULayer')
                        Li = ReluLayer.parse(L);
                    elseif isa(L, 'nnet.cnn.layer.MaxPooling2DLayer')
                        Li = MaxPooling2DLayer.parse(L);
                    elseif isa(L, 'nnet.cnn.layer.AveragePooling2DLayer')
                        Li = AveragePooling2DLayer.parse(L);
                    elseif isa(L, 'nnet.cnn.layer.FullyConnectedLayer')
                        Li = FullyConnectedLayer.parse(L);
                    else                     
                        fprintf('\nLayer %d is a %s which have not supported yet in nnv, please consider removing this layer for the analysis', i, class(L));
                        error('\nUnsupported Class of Layer');                     
                    end
                    
                    j = j + 1;
                    Ls{j} = Li;
                    
                end
                             
            end
            
            cnn = CNN(name, Ls, inputSize, outputSize);
            fprintf('\nParsing network is done successfully and %d Layers are neglected in the analysis phase', n - j);
            
        end
        
        
        
    end
    
    
    
end

