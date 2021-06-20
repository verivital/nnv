classdef MaxUnpooling2DLayer < handle
    % The MaxUnPooling 2D layer class in CNN
    %   Contain constructor and reachability analysis methods
    % Main references:
    % 1) 
    %    

    
    %   Dung Tran: 4/14/2020
    
    properties
        Name = 'max_unpooling_2d_layer';
        NumInputs = 3;
        InputNames = {'in',  'indices' , 'size'};
        NumOutputs = 1;
        OutputNames = {'out'};
        
        PairedMaxPoolingName = [];
    end
    
    
    methods
        
        % constructor of the class
        function obj = MaxUnpooling2DLayer(varargin)           
            % author: Dung Tran
            % date: 4/14/2020    
            % update: 
            
            switch nargin
                
                case 5
                    
                    name = varargin{1};
                    numInputs = varargin{2};
                    inputNames = varargin{3};
                    numOutputs = varargin{4};
                    outputNames = varargin{5};
                    
                    if ~ischar(name)
                        error('Name is not char');
                    end                    
                    
                    if numInputs < 1 
                        error('Invalid number of inputs');
                    end
                    
                    if ~iscell(inputNames)
                        error('InputNames should be a cell');
                    end
                    
                    if numOutputs < 1
                        error('Invalid number of outputs');
                    end
                    
                    if ~iscell(outputNames)
                        error('OutputNames should be a cell');
                    end
                    
                    
                    obj.Name = name;
                    obj.NumInputs = numInputs;
                    obj.InputNames = inputNames;
                    obj.NumOutputs = numOutputs;
                    obj.OutputNames = outputNames;
                    
                case 3
                    
                    name = varargin{1};
                    numInputs = varargin{2};
                    inputNames = varargin{3};
                    
                    if ~ischar(name)
                        error('Name is not char');
                    end                    
                    
                    if numInputs < 1 
                        error('Invalid number of inputs');
                    end
                    
                    if ~iscell(inputNames)
                        error('InputNames should be a cell');
                    end
                    
                    obj.Name = name;
                    obj.NumInputs = numInputs;
                    obj.InputNames = inputNames;
                    
                case 2
                    
                    name = 'max_unpooling_2d_layer';
                    numInputs = varargin{1};
                    inputNames = varargin{2};
                    
                    if ~ischar(name)
                        error('Name is not char');
                    end                    
                    
                    if numInputs < 1 
                        error('Invalid number of inputs');
                    end
                    
                    if ~iscell(inputNames)
                        error('InputNames should be a cell');
                    end
                    
                    obj.Name = name;
                    obj.NumInputs = numInputs;
                    obj.InputNames = inputNames;
                    
                case 0
                    
                    obj.Name = 'max_unpooling_2d_layer';
                                                        
                otherwise
                    error('Invalid number of inputs (should be 0 or 3)');
                                 
            end 
             
        end
        
        % set a paired max pooling name (used for analysis)
        function setPairedMaxPoolingName(obj, name)
            % @name: the name of the max pooling encoder
            
            % author: Dung Tran
            % date: 4/27/2020
            
            if ~ischar(name)
                error('Invalid name');
            else
                obj.PairedMaxPoolingName = name;
            end
            
        end
        
        
        % evaluation 
        function y = evaluate(~, input, indx, outputSize)
            % @input: input image
            % @indx: max index
            % @outputSize: inputSie
            % @y: output image
            
            % author: Dung Tran
            % date:4/19/2020
            
            dlX = dlarray(input, 'SSC');
            dlY = maxunpool(dlX, indx, outputSize); 
            y = extractdata(dlY);
            
        end
        
    end
    
    
    methods % reachability method
        
        function R = stepReachStar_singleInput(~, in_R, max_points, V, lb, ub)
            % @in_R: intermediate ImageStar
            % @max_points: max-point idexes
            % @V: the max-point star set: V = c + (a[1]*v[1] + a[2]*v[2] + ... + a[n]*v[n])
            % @lb: lower bound of the max point
            % @ub: upper bound of the max point
            
            % author: Dung Tran
            % date:4/28/2020
            
            
            N = size(max_points, 1); % number of local max points
            R = [];
            for l=1:N
                maxpoint = max_points(l,:,:);
                i = maxpoint(1);
                j = maxpoint(2);
                k = maxpoint(3);
                R1 = ImageStar(in_R.V, in_R.C, in_R.d, in_R.pred_lb, in_R.pred_ub, in_R.im_lb, in_R.im_ub);
                R1.MaxIdxs = in_R.MaxIdxs;
                R1.InputSizes = in_R.InputSizes;
                R1.V(i,j,k,:) = V;
                R1.im_lb(i, j, k) = lb;
                R1.im_ub(i, j, k) = ub;
                R = [R R1];
            end
        end
        
        
        
        function R = stepReachStar_multipleInputs(obj, in_R, max_points, V, lb, ub)
            % @in_R: an array of intermediate ImageStar
            % @max_points: max-point idexes
            % @V: max-point star set: V = c + (a[1]*v[1] + a[2]*v[2] + ... + a[n]*v[n])
            % @lb: lower bound of the max point
            % @ub: upper bound of the max point
            
            % author: Dung Tran
            % date:4/28/2020
            
            n = length(in_R);
            R = [];
            for i=1:n
                R = [R obj.stepReachStar_singleInput(in_R(i), max_points, V, lb, ub)];
            end
            
        end
        
        
        function R = reach_star(obj, IS)
            % @IS: input ImageStar
            % @R: output ImageStars
            
            % author: Dung Tran
            % date: 4/27/2020
            
            
            n = length(IS.MaxIdxs);
            newMaxIdxs = IS.MaxIdxs;
            newInputSizes = IS.InputSizes;
            for i=1:n
                if strcmp(IS.MaxIdxs{i}.Name, obj.PairedMaxPoolingName)
                    maxIdx = IS.MaxIdxs{i}.MaxIdx;
                    inputSize = IS.InputSizes{i}.InputSize;
                    newMaxIdxs{i} = [];
                    newInputSizes{i} = [];
                    break;
                end
            end
            
            numChannel = IS.numChannel;
            numPred = IS.numPred;
            h = IS.height;
            w = IS.width;
            if isempty(IS.im_lb)
                IS.estimateRanges;
            end
            
            % initial reach set
            V(:,:,numChannel, 1+ numPred) = zeros(inputSize(1), inputSize(2));
            im_lb(:,:,numChannel) = zeros(inputSize(1), inputSize(2));
            im_ub(:,:,numChannel) = zeros(inputSize(1), inputSize(2));
            R0 = ImageStar(V, IS.C, IS.d, IS.pred_lb, IS.pred_ub, im_lb, im_ub);
            R0.MaxIdxs = newMaxIdxs; 
            R0.InputSizes = newInputSizes;
            R = R0;
            for k=1:numChannel
                for i=1:h
                    for j=1:w
                        maxID = maxIdx{i, j, k}; % index of the local max points
                        n = size(maxID, 1);
                        maxPoints = [maxID k*ones(n,1)]; % the max point position: [height width channel]
                        V1 = IS.V(i,j,k,:);
                        lb = IS.im_lb(i, j, k);
                        ub = IS.im_ub(i ,j, k);
                        R = obj.stepReachStar_multipleInputs(R, maxPoints, V1, lb, ub);
                    end
                end
            end
            
        end
        
        
        % reach star with multiple inputs
        function IS = reach_star_multipleInputs(obj, in_images, option)
            % @in_images: an array of imagestar input sets
            % images: an array of imagestar output sets
            % option: '[]' or 'parallel'
            
            % author: Dung Tran
            % date: 4/28/2020
            
            
            n = length(in_images);
            IS = [];
            if strcmp(option, 'parallel')
                
                parfor i=1:n
                    IS = [IS obj.reach_star(in_images(i))];
                end
                
            elseif isempty(option) || strcmp(option, 'single')
                
                for i=1:n
                    IS = [IS obj.reach_star(in_images(i))];
                end
            else
                error('Unknown computation option');
            end
            
            
        end
        
                
    end
    
    methods % main reachability method 
        
        % general functio for reachability analysis
        function IS = reach(varargin)
            % @in_images: an input imagestar
            % @IS: output set
            % @option: = 'single' or 'parallel' 
            
            % author: Dung Tran
            % date: 4/28/2020
             
            switch nargin
                
                case 5
                    obj = varargin{1};
                    in_images = varargin{2};
                    method = varargin{3};
                    option = varargin{4};
                    %relaxFactor = varargin{5}; do not use
                
                case 4
                    obj = varargin{1};
                    in_images = varargin{2};
                    method = varargin{3};
                    option = varargin{4};
                
                case 3
                    obj = varargin{1};
                    in_images = varargin{2};
                    method = varargin{3};
                    option = [];
                
                otherwise
                    error('Invalid number of input arguments (should be 2, 3 or 4)');
            end
            
            if strcmp(method, 'approx-star')||strcmp(method, 'exact-star') || contains(method, 'relax-star')
                IS = obj.reach_star_multipleInputs(in_images, option);
            elseif strcmp(method, 'abs-dom')
                % abs-domain works similarly to approx-star method for unmax
                % pooling layer
                IS = obj.reach_star_multipleInputs(in_images, option);
            elseif strcmp(method, 'approx-zono')
                error('NNV havenot yet support approx-zono method');
            else
                error('Unknown reachability method');
            end
   
        end
        
        
    end
    
    
    methods(Static)
       
        
        
        % parse a trained MaxUnPooling2dLayer from matlab
        function L = parse(max_unpooling_2d_layer)
            % @max_unpooling_2d_Layer: an MaxUnPooling2DLayer from matlab deep
            % neural network tool box
            % @L : an MaxUnPooling2DLayer obj for reachability analysis purpose
            
            % author: Dung Tran
            % date: 4/14/2020
            
            
            if ~isa(max_unpooling_2d_layer, 'nnet.cnn.layer.MaxUnpooling2DLayer')
                error('Input is not a Matlab nnet.cnn.layer.MaxUnpooling2DLayer class');
            end
            
            L = MaxUnpooling2DLayer(max_unpooling_2d_layer.Name, max_unpooling_2d_layer.NumInputs, max_unpooling_2d_layer.InputNames, max_unpooling_2d_layer.NumOutputs, max_unpooling_2d_layer.OutputNames);
            fprintf('\nParsing a Matlab max pooling 2d layer is done successfully');
            
        end
        
    end
    
    
    
end

