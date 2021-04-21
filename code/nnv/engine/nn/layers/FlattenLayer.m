classdef FlattenLayer < handle
    % Flatten Layer object
    % Author: Dung Tran
    % Date: 6/9/2020
    
    properties
        Name = 'FlattenLayer';
        NumInputs = 1; % default
        InputNames = {'in'}; % default
        NumOutputs = 1; % default
        OutputNames = {'out'}; % default
        
        Type = ''; % type of flatten Layer, (from Matlab, or from Keras)
    end
    
    methods
        
        function obj = FlattenLayer(varargin)
            % @name: name of the layer
            % @NumInputs: number of inputs
            % @NumOutputs: number of outputs,
            % @InputNames: input names
            % @OutputNames: output names
            
            % author: Dung Tran
            % date:4/14/2020
            
            switch nargin
                
                case 5
                    name = varargin{1};
                    numInputs = varargin{2};
                    numOutputs = varargin{3};
                    inputNames = varargin{4};
                    outputNames = varargin{5};
                    
                case 1
                    name = varargin{1};
                    numInputs = 1;
                    numOutputs = 1;
                    inputNames = {'in'};
                    outputNames = {'out'};
                                        
                case 0
                    name = 'FlattenLayer';
                    numInputs = 1;
                    numOutputs = 1;
                    inputNames = {'in'};
                    outputNames = {'out'};
             
                otherwise
                    
                    error('Invalid number of input arguments, should be 0, 1, or 5');        
            end
            
            
            if ~ischar(name)
                error('Invalid name, should be a charracter array');
            end
            
            if numInputs < 1
                error('Invalid number of inputs');
            end
                       
            if numOutputs < 1
                error('Invalid number of outputs');
            end
            
            if ~iscell(inputNames)
                error('Invalid input names, should be a cell');
            end
            
            if ~iscell(outputNames)
                error('Invalid output names, should be a cell');
            end
            
            obj.Name = name;
            obj.NumInputs = numInputs;
            obj.NumOutputs = numOutputs;
            obj.InputNames = inputNames;
            obj.OutputNames = outputNames; 
                        
        end
            
    end
        
        
    methods
        
        % evaluate
        function flatten_im = evaluate(obj,image)
            %@image: an multi-channels image
            %@flatten_im: flatten image
            
            % author: Dung Tran
            % date: 6/9/2020
     
            n = size(image);
            
            if strcmp(obj.Type, 'nnet.keras.layer.FlattenCStyleLayer')
                
                if length(n) == 2
                    flatten_im = permute(image, [2 1]); % keras flatten layer
                    flatten_im = reshape(flatten_im, [1 n(1)*n(2)]);
                elseif length(n) == 3
                    flatten_im = permute(image, [3 2 1]); % keras flatten layer
                    flatten_im = reshape(flatten_im, [1 1 n(1)*n(2)*n(3)]);
                else
                    error('Invalid input image');
                end
                
            elseif strcmp(obj.Type, 'nnet.cnn.layer.FlattenLayer')        
                if length(n) == 2
                    flatten_im = reshape(image, [1 n(1)*n(2)]);
                elseif length(n) == 3
                    flatten_im = reshape(image, [1 1 n(1)*n(2)*n(3)]);
                else
                    error('Invalid input image');
                end      
            else
                error('Unknown type of flatten layer');
            end
            
                
        end
    end
    
    methods % reachability method
        
        function image = reach_single_input(obj, in_image)
            % @in_image: input imagestar
            % @image: output set
            
            % author: Dung Tran
            % date: 6/9/2020
            
            
            if ~isa(in_image, 'ImageStar') && ~isa(in_image, 'ImageZono')
                error('Input set is not an ImageStar or ImageZono');
            end
            
            N = in_image.height*in_image.width*in_image.numChannel;
                                   
            n = in_image.numPred;
            V(1, 1, :, in_image.numPred + 1) = zeros(N, 1);        
            for i=1:n+1
                V(1, 1,:,i) = obj.evaluate(in_image.V(:,:,:,i));
            end
            
            image = ImageStar(V, in_image.C, in_image.d, in_image.pred_lb, in_image.pred_ub);
            
        end
        
        % handle multiple inputs
        function S = reach_multipleInputs(obj, inputs, option)
            % @inputs: an array of ImageStars
            % @option: = 'parallel' or 'single'
            % @S: output ImageStar
            
            % author: Dung Tran
            % date: 1/6/2020
            
            n = length(inputs);
            if isa(inputs(1), 'ImageStar')
                S(n) = ImageStar;
            elseif isa(inputs(1), 'ImageZono')
                S(n) = ImageZono;
            else
                error('Unknown input data set');
            end
          
            if strcmp(option, 'parallel')
                parfor i=1:n
                    S(i) = obj.reach_single_input(inputs(i));
                end
            elseif strcmp(option, 'single') || isempty(option)
                for i=1:n
                    S(i) = obj.reach_single_input(inputs(i));
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
            % date: 6/9/2020
           
             
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
                    error('Invalid number of input arguments (should be 2, 3, 4, 5 or 6)');
            end
            
            if strcmp(method, 'approx-star') || strcmp(method, 'exact-star') || strcmp(method, 'abs-dom') || strcmp(method, 'approx-zono') || contains(method, "relax-star")
                IS = obj.reach_multipleInputs(in_images, option);
            else
                error('Unknown reachability method');
            end
  
        end
        
        
    end
    
    
    methods(Static)
        % parsing method
        
        function L = parse(layer)
            % @flatten_layer: 
            % @L: constructed layer
                        
            % author: Dung Tran
            % date: 4/14/2020
                      
            if ~isa(layer, 'nnet.keras.layer.FlattenCStyleLayer') && ~isa(layer, 'nnet.cnn.layer.FlattenLayer')
                error('Input is not a flatten layer');
            end
            
            L = FlattenLayer(layer.Name, layer.NumInputs, layer.NumOutputs, layer.InputNames, layer.OutputNames);
            fprintf('\nParsing a flatten layer is done successfully');
            L.Type = class(layer);
        end

    end
end

