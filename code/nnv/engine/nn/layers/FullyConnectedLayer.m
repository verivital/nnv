classdef FullyConnectedLayer < handle
    % The FullyConnectedLayer layer class in CNN
    %   Contain constructor and reachability analysis methods
    % Main references:
    % 1) An intuitive explanation of convolutional neural networks: 
    %    https://ujjwalkarn.me/2016/08/11/intuitive-explanation-convnets/
    % 2) More detail about mathematical background of CNN
    %    http://cs231n.github.io/convolutional-networks/
    %    http://cs231n.github.io/convolutional-networks/#pool
    % 3) Matlab implementation of Convolution2DLayer and MaxPooling (for training and evaluating purpose)
    %    https://www.mathworks.com/help/deeplearning/ug/layers-of-a-convolutional-neural-network.html
    %    https://www.mathworks.com/help/deeplearning/ref/nnet.cnn.layer.fullyconnectedlayer.html
    
    %   Dung Tran: 6/26/2019
    
    properties
        Name = 'fully_connected_layer';
        % Hyperparameters
        InputSize = 0;  % number of input
        OutputSize = 0; % number of output
        Weights = []; % weight matrix
        Bias  = []; % bias vector     
    end
    
    
    % constructor
    methods
        
        % constructor of the class
        function obj = FullyConnectedLayer(varargin)           
            % author: Dung Tran
            % date: 6/26/2019    
            % update: 
            
            switch nargin
                
                case 3
                    
                    name = varargin{1};
                    W = varargin{2};
                    b = varargin{3};
                                     
                    if ~ischar(name)
                        error('Name is not char');
                    else
                        obj.Name = name;
                    end
                    
                    if size(W, 1) ~= size(b, 1)
                        error('Inconsistent dimension between the weight matrix and bias vector');
                    end
                    
                    if size(b,2) ~= 1
                        error('Bias vector should have one column');
                    end
                    
                    obj.InputSize = size(W,2);
                    obj.OutputSize = size(W,1);
                    obj.Weights = W;
                    obj.Bias = b;
                    
                case 2
                    
                    W = varargin{1};
                    b = varargin{2};
                    
                    obj.Name = 'fully_connected_layer';
                    
                    if size(W, 1) ~= size(b, 1)
                        error('Inconsistent dimension between the weight matrix and bias vector');
                    end
                    
                    if size(b,2) ~= 1
                        error('Bias vector should have one column');
                    end
                    
                    obj.InputSize = size(W,2);
                    obj.OutputSize = size(W,1);
                    obj.Weights = W;
                    obj.Bias = b;
                           
                case 0
                    
                    obj.Name = 'fully_connected_layer';
                    % Hyperparameters
                    obj.InputSize = 0;
                    obj.OutputSize = 0; % step size for traversing input
                    obj.Weights = [];
                    obj.Bias = [];
                    
                otherwise
                    error('Invalid number of inputs (should be 0, 2 or 3)');
                                 
            end 
             
        end
        
    end
        
    % evaluation method
    methods
        
        function y = evaluate(obj, input)
            % @input: input
            % @y: output
            
            % author: Dung Tran
            % date: 6/26/2019
            
            
            n = size(input);
            
            N = 1;
            for i=1:length(n)
                N = N*n(i);
            end
            if N ~= obj.InputSize
                error('Inconsistency between the input dimension and InputSize of the network');
            end
            
            I = reshape(input, N, 1);
            y = double(obj.Weights)*I + double(obj.Bias);
            y = reshape(y, [1, 1, size(obj.Bias, 1)]);
             
        end 
       
    end   
     
    methods 
        
        %(reachability analysis using imagestar)
        function image = reach_star_single_input(obj, in_image)
            % @in_image: input imagestar
            % @image: output set
            
            % author: Dung Tran
            % date: 6/26/2019
            
            
            if ~isa(in_image, 'ImageStar')
                error('Input set is not an ImageStar');
            end
            
            N = in_image.height*in_image.width*in_image.numChannel;
            if N~= obj.InputSize
                error('Inconsistency between the size of the input image and the InputSize of the network');
            end
                       
            n = in_image.numPred;
            V(1, 1, :, in_image.numPred + 1) = zeros(obj.OutputSize, 1);        
            for i=1:n+1
                I = in_image.V(:,:,:,i);
                I = reshape(I,N,1); % flatten input
                if i==1
                    V(1, 1,:,i) = double(obj.Weights)*I + double(obj.Bias);
                else
                    V(1, 1,:,i) = double(obj.Weights)*I;
                end
            end
            
            image = ImageStar(V, in_image.C, in_image.d, in_image.pred_lb, in_image.pred_ub);
            
        end
        
        % handle multiple inputs
        function S = reach_star_multipleInputs(obj, inputs, option)
            % @inputs: an array of ImageStars
            % @option: = 'parallel' or 'single'
            % @S: output ImageStar
            
            % author: Dung Tran
            % date: 1/6/2020
            
            n = length(inputs);
            S(n) = ImageStar;
            if strcmp(option, 'parallel')
                parfor i=1:n
                    S(i) = obj.reach_star_single_input(inputs(i));
                end
            elseif strcmp(option, 'single') || isempty(option)
                for i=1:n
                    S(i) = obj.reach_star_single_input(inputs(i));
                end
            else
                error('Unknown computation option, should be parallel or single');
            end
            
        end
        
        
        
        %(reachability analysis using imagezono)
        function image = reach_zono(obj, in_image)
            % @in_image: input imagezono
            % @image: output set
            
            % author: Dung Tran
            % date: 1/2/2020
            
            
            if ~isa(in_image, 'ImageZono')
                error('Input set is not an ImageZono');
            end
            
            N = in_image.height*in_image.width*in_image.numChannels;
            if N~= obj.InputSize
                error('Inconsistency between the size of the input image and the InputSize of the network');
            end
                       
            n = in_image.numPreds;
            V(1, 1, :, in_image.numPreds + 1) = zeros(obj.OutputSize, 1);        
            for i=1:n+1
                I = in_image.V(:,:,:,i);
                I = reshape(I,N,1);
                if i==1
                    V(1, 1,:,i) = double(obj.Weights)*I + double(obj.Bias);
                else
                    V(1, 1,:,i) = double(obj.Weights)*I;
                end
            end
            
            image = ImageZono(V);
            
        end
        
        % handle mulitples inputs
        function S = reach_zono_multipleInputs(obj, inputs, option)
            % @inputs: an array of ImageZonos
            % @option: = 'parallel' or 'single'
            % @S: output ImageZono
            
            % author: Dung Tran
            % date: 1/6/2020
            
            n = length(inputs);
            S(n) = ImageZono;
            if strcmp(option, 'parallel')
                parfor i=1:n
                    S(i) = obj.reach_zono(inputs(i));
                end
            elseif strcmp(option, 'single') || isempty(option)
                for i=1:n
                    S(i) = obj.reach_zono(inputs(i));
                end
            else
                error('Unknown computation option, should be parallel or single');
            end
            
        end
        
        
        
        % reachability analysis with multiple inputs
        function IS = reach(varargin)
            % @in_image: an input imagestar
            % @image: output set
            % @option: = 'single' or 'parallel' 
            
            % author: Dung Tran
            % date: 6/26/2019
            % update: 1/6/2020  update reason: add zonotope method
             
            switch nargin
                
                 case 7
                    obj = varargin{1};
                    in_images = varargin{2};
                    method = varargin{3};
                    option = varargin{4};
                    % relaxFactor = varargin{5}; do not use
                    % dis_opt = varargin{6}; do not use
                    % lp_solver = varargin{7}; do not use
                
                case 6
                    obj = varargin{1};
                    in_images = varargin{2};
                    method = varargin{3};
                    option = varargin{4};
                    %relaxFactor = varargin{5}; do not use
                    % dis_opt = varargin{6}; do not use
                
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
                    option = varargin{4}; % computation option

                case 3
                    obj = varargin{1};
                    in_images = varargin{2}; % don't care the rest inputs
                    method = varargin{3};
                    option = [];
                otherwise
                    error('Invalid number of input arguments (should be 2, 3, 4, 5, or 6)');
            end
            
            if strcmp(method, 'approx-star') || strcmp(method, 'exact-star') || strcmp(method, 'abs-dom') || contains(method, "relax-star")
                IS = obj.reach_star_multipleInputs(in_images, option);
            elseif strcmp(method, 'approx-zono')
                IS = obj.reach_zono_multipleInputs(in_images, option);
            else
                error('Unknown reachability method');
            end
            
            
            
        end
        
    end
    
    
    methods(Static)
        
         
        % parse a trained averagePooling2dLayer from matlab
        function L = parse(fully_connected_layer)
            % @fully_connecteted_Layer: a fully connected layer from matlab deep
            % neural network tool box
            % @L : a FullyConnectedLayer obj for reachability analysis purpose
            
            % author: Dung Tran
            % date: 6/26/2019
            
            
            if ~isa(fully_connected_layer, 'nnet.cnn.layer.FullyConnectedLayer')
                error('Input is not a Matlab nnet.cnn.layer.FullyConnectedLayer class');
            end
                        
            L = FullyConnectedLayer(fully_connected_layer.Name, fully_connected_layer.Weights, fully_connected_layer.Bias);         
            fprintf('\nParsing a Matlab fully connected layer is done successfully');
            
        end
        
        
    end
    
    
    
    
end

