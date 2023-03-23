classdef DepthConcatenationLayer < handle
    % Depth Concatenation Layer object
    % Concatenate arrays along 3rd dimension (channels)
    % Author: Diego Manzanas Lopez
    % Date: 03/07/2023
    
    properties
        Name = 'DepthConcatLayer';
        NumInputs = 1; % default
        InputNames = {'in'}; % default
        NumOutputs = 1; % default
        OutputNames = {'out'}; % default
    end
    
    methods
        
        function obj = DepthConcatenationLayer(varargin)
            % @name: name of the layer
            % @NumInputs: number of inputs
            % @NumOutputs: number of outputs,
            % @InputNames: input names
            % @OutputNames: output names
            
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
                    name = 'DepthConcatLayer';
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
        
        
    methods % main methods
        
        % evaluate (TODO)
        function outputs = evaluate(obj, inputs)
            % depth_concatenation layers takes usually two inputs, but allows many (N)
            % first input is obj, the rest are arrays from different layers
            % Initialize image as the first one
            outputs = inputs{1};
            % Concatenate the inputs 
            for k=2:length(inputs)
                outputs = cat(obj.Dim, outputs, inputs{k});
            end
                
        end
 
        % reach (TODO)
        function outputs = reach_single_input(obj, inputs)
            % @inputs: input imagestar from each connected layer
            % @outputs: output set
            
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
                      
            if ~isa(layer, 'nnet.cnn.layer.DepthConcatenationLayer')
                error('Input is not a depth concatenation layer');
            end
            
            L = DepthConcatenationLayer(layer.Name, layer.NumInputs, layer.NumOutputs, layer.InputNames, layer.OutputNames);
            fprintf('\nParsing a depth concatenation layer is done successfully');
        end

    end
end

