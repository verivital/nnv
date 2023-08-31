classdef GlobalAveragePooling1DLayer < handle
    % GlobalAveragePooling 2D Layer object
    % downsamples the input by computing the mean of the height and width dimensions
    % https://www.mathworks.com/help/deeplearning/ref/nnet.cnn.layer.globalAveragePooling2dLayer.html
    % 
    % Author: Neelanjana Pal
    % Date: 06/07/2023
    
    properties
        Name = 'GlobalAvgPooling1DLayer';          % default
        NumInputs = 1;              % default
        InputNames = {'in1'}; % default
        NumOutputs = 1;             % default
        OutputNames = {'out'};      % default
    end
    
    methods % constructor
        
        % create layer
        function obj = GlobalAveragePooling1DLayer(varargin)
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
        function output = evaluate(obj, input)
            % addition layer takes usually two inputs, but allows many (N)
            %
            input = dlarray(input',"SC");
            output = extractdata(avgpool(input,'global'))';
        end
        
        function output = evaluateSequence(obj, input)
            output = obj.evaluate(input);
        end

        % reach (TODO)
        function image = reach_single_input(obj, input)
            % @in_image: input imagestar
            % @image: output set
            
            if ~isa(input, 'ImageStar')
                error('The input is not an ImageStar object');
            end
            % Y = obj.evaluate(input.V);                       
            % output = ImageStar(Y, input.C, input.d, input.pred_lb, input.pred_ub);
            if isempty(input.im_lb) && isempty(input.im_ub)
                c = obj.evaluateSequence(input.V(:,:,:,1));
                layer = obj;
                parfor i = 2: size(in_image.V,4)
                    V(:,:,:,i) = layer.evaluateSequence(input.V(:,:,:,i));
                end
                V(:,:,:,1) = c;
                image = ImageStar(V, input.C, input.d, input.pred_lb, input.pred_ub);
            else
                im_lb = input.im_lb;
                im_ub = input.im_ub;
                lb = obj.evaluateSequence(im_lb);
                ub = obj.evaluateSequence(im_ub);
    
                image = ImageStar(lb,ub);
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
                IS = obj.reach_multipleInputs(in_images, option);
            else
                error('Unknown reachability method');
            end
  
        end
        
        function IS = reachSequence(varargin)
            obj = varargin{1};
            IS = obj.reach(varargin{2:end});
        end
    end
    
    
    methods(Static)
        
        % parsing method
        function L = parse(layer)
            % create NNV layer from matlab
                      
            if ~isa(layer, 'nnet.cnn.layer.GlobalAveragePooling1DLayer') 
                error('Input is not a GlobalAveragePooling1DLayer layer');
            end
            L = GlobalAveragePooling1DLayer(layer.Name, layer.NumInputs, layer.NumOutputs, layer.InputNames, layer.OutputNames);
        end

    end
end

