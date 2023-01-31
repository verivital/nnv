classdef ReshapeLayer < handle
    %ReshapeLayer in NN
    % reshapes the input into a specific dimensions, tipically a 2D into
    % 4D, from FullyConnected to Conv2D layers

    % author: Diego Manzanas
    % Date: 12/08/2022
    
    properties
        Name = 'reshape_layer'; % default
        NumInputs = 1;          % default
        InputNames = {'in'};    % default
        NumOutputs = 1;         % default
        OutputNames = {'out'};  % default
        targetDim = [];         % default
    end
    
    methods
        function obj = ReshapeLayer(varargin)
            switch nargin
                
                case 6
                    name = varargin{1};
                    numInputs = varargin{2};
                    numOutputs = varargin{3};
                    inputNames = varargin{4};
                    outputNames = varargin{5};
                    targetDim = varargin{6};
                    
                case 2
                    name = varargin{1};
                    targetDim = varargin{2};
                    numInputs = 1;
                    numOutputs = 1;
                    inputNames = {'in'};
                    outputNames = {'out'};
                                        
                case 0
                    name = 'reshape_layer';
                    numInputs = 1;
                    numOutputs = 1;
                    inputNames = {'in'};
                    outputNames = {'out'};
                    targetDim = [];
             
                otherwise
                    
                    error('Invalid number of input arguments, should be 0, 2, or 6');        
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
            obj.targetDim = targetDim;
                        
        end
    end

    methods
        
        % evaluate
        function reshape_x = evaluate(obj,image)
            %@image: an multi-channels image
            %@flatten_im: flatten image
     
            reshape_x = reshape(image, obj.targetDim);
                
        end
    end
    
    methods % reachability method
        
        function image = reach_single_input(obj, in_image)
            % @in_image: input imagestar
            % @image: output set
            
            % TODO: implement this function, just need to modify the
            % dimensions of an ImageStar or convert a Star to an ImageStar
            % Should also support ImageZono and Zono
            error("TODO, Working on adding support for this layer.")
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
            % @layer: Custom Reshape Layer (created during importONNXLayers)
            % @L: constructed layer
                      
            if ~contains(class(layer), 'ReshapeLayer')
                error('Input is not a reshape layer');
            end
            
            params = layer.ONNXParams.Nonlearnables;
            par_fields = fields(params);
            if length(par_fields) == 1
                params = struct2cell(params);
                targetDim = extractdata(params{1});
                targetDim = reshape(targetDim, [1 length(targetDim)]);
            else
                error('Parsing Reshape Layer was unsuccessful. We only support reshape layer with one Nonlearnable parameter.')
            end

            L = ReshapeLayer(layer.Name, layer.NumInputs, layer.NumOutputs, layer.InputNames, layer.OutputNames, targetDim);
            fprintf('Parsing a reshape layer is done successfully \n');
        end

    end

end

