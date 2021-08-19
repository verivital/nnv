classdef ElementwiseAffineLayer < handle
    % The ElementwiseAffineLayer layer class in CNN
    % author: Neelanjana Pal
    % date: 6/28/2021
    
    properties
        Name = 'elementwise_affine_layer';
        % Hyperparameters
        NumInputs = 1; % default
        InputNames = {'in'}; % default
        NumOutputs = 1; % default
        OutputNames = {'out'}; % default
        
        Scale = []; % Scale vector
        Offset  = []; % Offset vector
    end
    
    
    % constructor
    methods
        
        % constructor of the class
        function obj = ElementwiseAffineLayer(varargin)           
            % author: Neelanjana Pal
            % date: 6/28/2021     
            
            switch nargin
                
                case 7
                    obj.Name = varargin{1};
                    obj.NumInputs = varargin{2};
                    obj.InputNames = varargin{3};
                    obj.NumOutputs = varargin{4};
                    obj.OutputNames = varargin{5};
                    obj.Scale = varargin{6};
                    obj.Offset = varargin{7};
                    %prevL = varargin{4};
                case 3
                    obj.Name = varargin{1} ; 
                    obj.Scale = varargin{2};
                    obj.Offset = varargin{3};
                case 2
                    obj.Scale = varargin{1};
                    obj.Offset = varargin{2};
                    
                    obj.Name = 'elementwise_affine_layer';              
                    
                case 0
                    
                    obj.Name = 'elementwise_affine_layer';
                           
                otherwise
                    error('Invalid number of inputs (should be 0 or 1)');
            end 
             
        end
        
    end
        
    % evaluation method
    methods
        
        function y = evaluate(obj, input)
            % @input: 2 or 3-dimensional array, for example, input(:, :, :), 
            % @y: 2 or 3-dimensional array, for example, y(:, :, :
            
            % author: Neelanjana Pal
            % date: 6/28/2021
            
            y = input;
            % assuming Scale or Object is with dim 1x1xinput.numChannel
            if ~isscalar(obj.Scale) && size(obj.Scale, 3) ~= size(input,3)%in_image.numChannel
                error('Inconsistent number of channels between Scale array and the ImageStar');
            elseif ~isscalar(obj.Scale) && size(obj.Scale, 3) == size(input,3)%in_image.numChannel
                for i=1:size(obj.Scale, 3)
                    y(:,:,i)=input(:,:,i)*obj.Scale(i);
                end
            end
            
            if ~isscalar(obj.Offset) && size(obj.Offset, 3) ~= size(input,3)%in_image.numChannel
                error('Inconsistent number of channels between Offset array and the ImageStar');
            elseif ~isscalar(obj.Offset) && size(obj.Offset, 3) == size(input,3)%in_image.numChannel
                for i=1:size(obj.Offset, 3)
                    y(:,:,i)=input(:,:,i)+obj.Offset(i);
                end
            end
        end 
       
    end   
     
 methods % reachability method
        
        function image = reach_single_input(obj, in_image)
            % @in_image: input imagestar
            % @image: output set
            
            % author: Neelanjana Pal
            % date: 6/28/2021
            
            
            if ~isa(in_image, 'ImageStar') && ~isa(in_image, 'ImageZono')
                error('Input set is not an ImageStar or ImageZono');
            end
            
            if ~isempty(obj.Scale) || ~isempty(obj.Offset)
                new_V = obj.evaluate(in_image.V);
            else
                new_V = in_image.V;
            end
                     
            image = ImageStar(new_V, in_image.C, in_image.d, in_image.pred_lb, in_image.pred_ub);            
            
        end
        
        % handle multiple inputs
        function S = reach_multipleInputs(obj, inputs, option)
            % @inputs: an array of ImageStars
            % @option: = 'parallel' or 'single'
            % @S: output ImageStar
            
            % author: Neelanjana Pal
            % date: 6/28/2021
            
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
            
            % author: Neelanjana Pal
            % date: 6/28/2021
           
             
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
        
         
        % parse a trained elementwise affine layer from matlab
        function L = parse(elementwise_affine_layer)
            % @elementwise_affine_layer: a elementwise affine layer from matlab deep
            % neural network tool box
            % @L : a ElementwiseAffineLayer obj for reachability analysis purpose
            
            % author: Neelanjana Pal
            % date: 6/28/2021
            
            
            if ~isa(elementwise_affine_layer, 'nnet.onnx.layer.ElementwiseAffineLayer')
                error('Input is not a Matlab nnet.onnx.layer.ElementwiseAffineLayer class');
            end
                        
            L = ElementwiseAffineLayer(elementwise_affine_layer.Name, elementwise_affine_layer.Scale, elementwise_affine_layer.Offset);         
            fprintf('\nParsing a Matlab elementwise affine layer is done successfully');
            
        end
        
        
    end
    
    
    
    
end

