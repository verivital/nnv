classdef ConcatenationLayer < handle
    % Concatenation Layer object
    % Concatenate arrays along a specific dimension
    % Typically use to concatenate 1D or 2D features  (vector/matrices, rather than images such as in DepthConcatenation)
    % see: https://www.mathworks.com/help/deeplearning/ref/nnet.cnn.layer.concatenationlayer.html
    % 
    % Author: Diego Manzanas Lopez
    % Date: 03/15/2023
    
    properties
        Name = 'ConcatLayer';
        NumInputs = 1; % default
        InputNames = {'in'}; % default
        NumOutputs = 1; % default
        OutputNames = {'out'}; % default
        Dim = 0 % unspecified, must be a positive integer
    end
    
    methods % constructor
        
        % create layer
        function obj = ConcatenationLayer(varargin)
            % @name: name of the layer
            % @NumInputs: number of inputs
            % @NumOutputs: number of outputs,
            % @InputNames: input names
            % @OutputNames: output names
            % @Dim: dimension for concatenating inputs
            
            switch nargin
                
                case 6
                    name = varargin{1};
                    numInputs = varargin{2};
                    numOutputs = varargin{3};
                    inputNames = varargin{4};
                    outputNames = varargin{5};
                    dim = varargin{6};
                case 2
                    name = varargin{1};
                    dim = varargin{2};
                    numInputs = 1;
                    numOutputs = 1;
                    inputNames = {'in'};
                    outputNames = {'out'};
                otherwise
                    error('Invalid number of input arguments, should be 2 or 6');        
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

            if ~iscalar(dim) || dim < 1
                error('Concat dimension must be a positive integer')
            end
            
            obj.Name = name;
            obj.NumInputs = numInputs;
            obj.NumOutputs = numOutputs;
            obj.InputNames = inputNames;
            obj.OutputNames = outputNames; 
            obj.Dim = dim;
        end
            
    end
        
        
    methods % main methods
        
        % evaluate
        function outputs = evaluate(obj, inputs)
            % concatenation layers takes usually two inputs, but allows many (N)
            %
            % first input is obj, the second is a cell array containing the inputs to concatenate
            % Initialize outputs as the first one
            outputs = inputs{1};
            % Concatenate the inputs 
            for k=2:length(inputs)
                outputs = cat(obj.Dim, outputs, inputs{k});
            end
                
        end
 
        % reach (TODO)
        function image = reach_single_input(obj, in_image)
            % @in_image: input imagestar
            % @image: output set
            
        end
        
        % handle multiple inputs
        function S = reach_multipleInputs(obj, inputs, option)
            % @inputs: an array of ImageStars
            % @option: = 'parallel' or 'single'
            % @S: output ImageStar
            
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
            % create NNV layer from matlab
                      
            if ~isa(layer, 'nnet.cnn.layer.ConcatenationLayer')
                error('Input is not a concatenation layer');
            end
            L = ConcatenationLayer(layer.Name, layer.NumInputs, layer.NumOutputs, layer.InputNames, layer.OutputNames, layer.Dimension);
            fprintf('\nParsing a concatenation layer is done successfully');
        end

    end
end

