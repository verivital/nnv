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
                        if ~isa(L, 'ImageInputLayer') && ~isa(L, 'BatchNormalizationLayer') && ~isa(L, 'AveragePooling2DLayer') && ~isa(L, 'Conv2DLayer') && ~isa(L, 'FullyConnectedLayer') && ~isa(L, 'MaxPooling2DLayer') && ~isa(L, 'ReluLayer')
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
        function y = evaluate(obj, x)
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
    
    
    methods % clasification and verification methods
        
        function label_id = classify(varargin)
            % @in_image: an image or an imagestar
            % @label_id: output index of classified object
            
            % author: Dung Tran
            % date: 8/21/2019
            
            switch nargin
                case 2
                    obj = varargin{1};
                    in_image = varargin{2};
                    
                case 3
                    obj = varargin{1};
                    in_image = varargin{2};
                    method = varargin{3};
                    numOfCores = 1;
                    
                case 4
                    
                    obj = varargin{1};
                    in_image = varargin{2};
                    method = varargin{3};
                    numOfCores = varargin{4};
                  
                otherwise
                    error('Invalid number of inputs, should be 1, 2, or 3');
            end
            
                     
            
            if ~isa(in_image, 'ImageStar')
                [y, ~] = obj.evaluate(in_image);
                y = reshape(y, [obj.OutputSize, 1]);
                [~, label_id] = max(y); 
            else

                if nargin == 2
                    method = 'approx-star'; % default reachability method
                    numOfCores = 1; 
                end
                fprintf('\n=============================================');
                obj.reach(in_image, method, numOfCores);

                RS = obj.reachSet{obj.numLayers + 1};
                n = length(RS);
                label_id = cell(n, 1);
                for i=1:n
                    rs = RS(i);                        
                    new_rs  = ImageStar.reshape(rs, [obj.OutputSize 1 1]);                    
                    max_id = new_rs.get_localMax_index([1 1], [obj.OutputSize 1], 1);
                    label_id{i} = max_id(:, 1);
                end

            end
 
        end
        
        
        function [robust, counterExamples] = verifyRobustness(varargin)
            % @robust: = 1: the network is robust
            %          = 0: the network is notrobust
            %          = 2: robustness is uncertain
            % @counterExamples: a set of counter examples 
            % 
            % author: Dung Tran
            % date: 8/21/2019
            
            switch nargin
                
                case 3
                    obj = varargin{1};
                    in_image = varargin{2};
                    correct_id = varargin{3};
                    method = 'approx-star';
                    numOfCores = 1;
                    
                case 4
                    obj = varargin{1};
                    in_image = varargin{2};
                    correct_id = varargin{3};
                    method = varargin{4};
                    numOfCores = 1; 
                    
                case 5
                    obj = varargin{1};
                    in_image = varargin{2};
                    correct_id = varargin{3};
                    method = varargin{4};
                    numOfCores = varargin{5}; 
                    
                otherwise
                    error('Invalid number of inputs, should be 2, 3, or 4');
                     
            end
                        
            label_id = obj.classify(in_image, method, numOfCores);      
            n = length(label_id); 
            % check the correctness of classifed label
            counterExamples = [];
            incorrect_id_list = []; 
            for i=1:n
                ids = label_id{i};
                m = length(ids);
                id1 = [];
                for j=1:m
                    if ids(j) ~= correct_id
                       id1 = [id1 ids(j)];
                    end
                end
                incorrect_id_list = [incorrect_id_list id1];             
                
                % construct counter example set
                if ~isempty(id1)
                    
                    if ~isa(in_image, 'ImageStar')
                        
                        counterExamples = in_image; 
                    
                    elseif strcmp(method, 'exact-star')
                        
                        rs = obj.reachSet{obj.numLayers + 1}(i);
                        L = length(id1); 
                        for l=1:L                    
                            [new_C, new_d] = ImageStar.addConstraint(rs, [1, 1, correct_id], [1, 1, id1(l)]);  
                            counter_IS = ImageStar(in_image.V, new_C, new_d, in_image.pred_lb, in_image.pred_ub);
                            counterExamples = [counterExamples counter_IS];
                        end
                        
                    end
   
                end
      
            end
            
            
            if isempty(incorrect_id_list)
                robust = 1;
                fprintf('\n=============================================');
                fprintf('\nTHE NETWORK IS ROBUST');
                fprintf('\nClassified index: %d', correct_id);
                
            else
                
                if strcmp(method, 'exact-star')
                    robust = 0;
                    fprintf('\n=============================================');
                    fprintf('\nTHE NETWORK IS NOT ROBUST');
                    fprintf('\nLabel index: %d', correct_id);
                    fprintf('\nClassified index:');
                    
                    n = length(incorrect_id_list);
                    for i=1:n
                        fprintf('%d ', incorrect_id_list(i));
                    end
                    
                else
                    robust = 2;
                    fprintf('\n=============================================');
                    fprintf('\nThe robustness of the network is UNCERTAIN due to the conservativeness of approximate analysis');
                    fprintf('\nPossible classified index: ');
                                       
                    n = length(incorrect_id_list);
                    for i=1:n
                        fprintf('%d ', incorrect_id_list(i));
                    end                   
                    fprintf('\nPlease try to verify the robustness with exact-star (exact analysis) option');
                end
                
            end
            

        end
        
        
    end
    
    methods % evaluate robustness
        
        % evaluate robustness of a network on an array of input (test) sets
        function r = evaluateRobustness(varargin)
            % @in_images: input sets
            % @correct_ids: an array of correct labels corresponding to the input sets
            % @method: reachability method: 'exact-star', 'approx-star',
            % 'approx-zono' and 'abs-dom'
            % @numCores: number of cores used in computation
            % @r: robustness value (in percentage)           
            
            % author: Dung Tran
            % date:1/9/2020
            
            switch nargin
                case 5
                    obj = varargin{1};
                    in_images = varargin{2};
                    correct_ids = varargin{3};
                    method = varargin{4};
                    numOfCores = varargin{5};
                case 4
                    obj = varargin{1};
                    in_images = varargin{2};
                    correct_ids = varargin{3};
                    method = varargin{4};
                    numOfCores = 1;
                otherwise
                    error('Invalid number of input arguments, should be 3 or 4');                    
            end
            
            
            N = length(in_images);
            if length(correct_ids) ~= N
                error('Inconsistency between the number of correct_ids and the number of input sets');
            end
            
                      
            count = zeros(1, N);
            if ~strcmp(method, 'exact-star')                
                % compute reachable set
                outputSets = obj.reach(in_images, method, numOfCores);
                
                % verify reachable set              
                if numOfCores> 1
                    parfor i=1:N
                        count(i) = CNN.isRobust(outputSets(i), correct_ids(i));                            
                    end
                else
                    for i=1:N
                        count(i) = CNN.isRobust(outputSets(i), correct_ids(i));
                    end
                end
                
            end
            
            if strcmp(method, 'exact-star')
                % compute reachable set
                if numOfCores > 1
                    parfor i=1:N
                        outputSets = obj.reach(in_images(i), method);
                        % verify reachable set
                        M = length(outputSets);
                        count1 = 0;
                        for j=1:M
                            count1 = count1 + CNN.isRobust(outputSets(j), correct_ids(i));
                        end
                        if count1 == M
                            count(i) = 1;
                        end
                        
                    end
                else
                    for i=1:N
                        outputSets = obj.reach(in_images(i), method);
                        % verify reachable set
                        M = length(outputSets);
                        count1 = 0;
                        for j=1:M
                            count1 = count1 + CNN.isRobust(outputSets(j), correct_ids(i));
                        end
                        if count1 == M
                            count(i) = 1;
                        end
                    end
                end
                                    
            end          
                        
            r = sum(count)/N; 
            
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
                    elseif isa(L, 'nnet.cnn.layer.BatchNormalizationLayer')
                        Li = BatchNormalizationLayer.parse(L);
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
        
             
        % check robustness using the outputSet
        function bool = isRobust(outputSet, correct_id)
            % @outputSet: the outputSet we need to check
            % @correct_id: the correct_id of the classified output
            
            % author: Dung Tran
            % date: 1/11/2020
            
            if correct_id > outputSet.numChannel || correct_id < 1
                error('Invalid correct id');
            end
            
            count = 0;
            for i=1:outputSet.numChannel
                if correct_id ~= i
                    if outputSet.is_p1_larger_p2([1 1 i], [1 1 correct_id]);
                       bool = 0;
                       break;
                    else
                        count = count + 1;
                    end
                end
            end
            
            if count == outputSet.numChannel - 1
                bool = 1;
            end
            
        end
        
        
        % classify outputset
        function classified_id = classifyOutputSet(outputSet)
            % @outputSet: is the output of a CNN, it is an ImageStar or
            %             ImageZono object
            % @classified_id: the classified_id = the id of classified
            % output or [] if we cannot classified outputset
            % the classified id is corresponding to the output that has
            % maximum value
            
            
            % author: Dung Tran
            % date: 1/10/2020
            
            if ~isa(outputSet, 'ImageStar') && ~isa(outputSet, 'ImageZono')
                error('Output set is not an ImageStar or an ImageZono');
            end           
  
            [lb, ub] = outputSet.getRanges; 
     
            
            [max_lb, max_lb_id] = max(lb);
            ub(max_lb_id) = [];
            ub = (ub > max_lb); 
            
            if sum(ub) == 0
                classified_id = max_lb_id;
            else
                classified_id = [];
            end
            
        end
        
        
        
        
    end
    
    
    
end

