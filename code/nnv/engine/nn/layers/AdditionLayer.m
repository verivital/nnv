classdef AdditionLayer < handle
    % Addition Layer object
    % Adds inputs from multiple neural network layers element-wise
    % https://www.mathworks.com/help/deeplearning/ref/nnet.cnn.layer.additionlayer.html
    % 
    % Author: Neelanjana Pal
    % Date: 06/07/2023
    
    properties
        Name = 'AddLayer';          % default
        NumInputs = 1;              % default
        InputNames = {'in1'}; % default
        NumOutputs = 1;             % default
        OutputNames = {'out'};      % default
    end
    
    methods % constructor
        
        % create layer
        function obj = AdditionLayer(varargin)
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
                    inputNames = {'in1'};
                    outputNames = {'out'};
                otherwise
                    error('Invalid number of input arguments, should be 1 or 5');        
            end
            
            if ~ischar(name)
                error('Invalid name, should be a charracter array');
            end
            
            if numInputs < 2
                error('Invalid number of inputs');
            end
                       
            if numOutputs < 1
                error('Invalid number of outputs');
            end
            
            if ~iscell(inputNames)
                error('Invalid input names, should be a cell');
            end
            
            if ~iscell(outputNames)
                error('num of inputs do not match with num of input names');
            end
            
            obj.Name = name;
            obj.NumInputs = numInputs;
            obj.NumOutputs = numOutputs;
            obj.InputNames = inputNames;
            obj.OutputNames = outputNames; 
        end
            
    end
        
        
    methods % main methods
        
        % evaluate
        function outputs = evaluate(obj, inputs)
            % addition layer takes usually two inputs, but allows many (N)
            %
            % first input is obj, the second is a cell array containing the
            % inputs to add
            % Initialize outputs as the first one
            outputs = inputs{1};
            % add the inputs 
            for k=2:length(inputs)
                outputs = outputs + inputs{k};
            end
                
        end
 
        % reach (TODO)
        function outputs = reach_single_input(obj, inputs)
            % @in_image: input imagestar
            % @image: output set
            
            outputs = inputs{1};
            for k = 2 : length(inputs)
                outputs = outputs.MinkowskiSum(inputs{k});
            end
            
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
%                 IS = obj.reach_multipleInputs(in_images, option);
                IS = obj.reach_single_input(in_images);
            else
                error('Unknown reachability method');
            end
  
        end
        
    end
    
    
    methods(Static)
        
        % parsing method
        function L = parse(layer)
            % create NNV layer from matlab
                      
            if ~isa(layer, 'nnet.cnn.layer.AdditionLayer') 
                error('Input is not a AdditionLayer layer');
            end
            L = AdditionLayer(layer.Name, layer.NumInputs, layer.NumOutputs, layer.InputNames, layer.OutputNames);
%             fprintf('\nParsing a addition layer is done successfully');
        end

    end
end

