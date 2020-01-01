classdef BatchNormalizationLayer < handle
    % The Batch Normalization Layer class in CNN
    %   Contain constructor and reachability analysis methods   
    %   Dung Tran: 1/1/2020
    
    properties
        Name = 'BatchNormalizationLayer';
        
        TrainedMean = [];
        TrainedVariance = [];
        
        % Hyperparameters 
        Epsilon = 0.00001;  % default value
        
        % Learnable parameters
        Offset = [];
        Scale = [];
               
    end
    
    
    % setting hyperparameters method
    methods
        
        % constructor of the class
        function obj = BatchNormalizationLayer(varargin)           
            % author: Dung Tran
            % date: 1/1/2020    
            % update: 
            
            
            if mod(nargin, 2) ~= 0
                error('Invalid number of arguments');
            end
            
            for i=1:nargin-1
                
                if mod(i, 2) == 0
                    
                    if strcmp(varargin{i}, 'Name')
                        obj.Name = varargin{i+1};
                    elseif strcmp(varargin{i}, 'TrainedMean')
                        obj.TrainedMean = varargin{i+1};
                    elseif strcmp(varargin{i}, 'TrainedVariance')
                        obj.TrainedVariance = varargin{i+1};
                    elseif strcmp(varargin{i}, 'Epsilon')
                        obj.Epsilon = varargin{i+1};
                    elseif strcmp(varargin{i}, 'Offset')
                        obj.Offset = varargin{i+1};
                    elseif strcmp(varargin{i}, Scale)
                        obj.Scale = varargin{i+1};
                    end
                    
                end
                
            end
                
             
        end
        
        
    end
        
    % evaluation method
    methods
        
        function y = evaluate(obj, input)
            % @input: input image
            % @y: output image with normalization
            
            % author: Dung Tran
            % date: 1/1/2020
                             
            if ~isempty(obj.TrainedMean) && ~isempty(obj.TrainedVariance) && ~isempty(obj.Epsilon) && ~isempty(obj.Offset) && ~isempty(obj.Scale)
                x1 = double(input) - obj.TrainedMean;
                x1 = x1/(sqrt(obj.TrainedVariance^2 + obj.Epsilon));
                y = obj.Scale*x1 + obj.Offset;
            else
                y = double(input);
            end
                               
        end
        
        
    end
        
    % exact reachability analysis using star set
    methods
        
        function images = reach(varargin)
            % @in_image: an input imagestar
            % @image: output set
            % @option: = 'single' or 'parallel' 
            
            % author: Dung Tran
            % date: 6/26/2019
             
            switch nargin
                
                case 4
                    obj = varargin{1};
                    in_images = varargin{2};
                    method = varargin{3};
                    option = varargin{4};
                
                case 3
                    obj = varargin{1};
                    in_images = varargin{2};
                    option = varargin{3};
                case 2
                    obj = varargin{1};
                    in_images = varargin{2};
                    option = [];
                otherwise
                    error('Invalid number of input arguments (should be 1, 2 or 3)');
            end
            
            
            n = length(in_images);
            for i=1:n
                if ~isa(in_images(i), 'ImageStar')
                    error('The %d^th input is not an ImageStar', i);
                end
            end
            
            images(n) = ImageStar;
            if strcmp(option, 'parallel')
                parfor i=1:n
                    V = double(in_images(i).V); % convert to double precision
                    x = in_images(i).affineMap(1/sqrt(obj.TrainedVariance^2 + obj.Epsilon), []);
                    x = x.affineMap(1, obj.TrainedMean);
                    y = 
                    
                    
                    temp = obj.evaluate(V(:,:,:,1));
                    V(:,:,:,1) = temp;
                    images(i) = ImageStar(V, in_images(i).C, in_images(i).d, in_images(i).pred_lb, in_images(i).pred_ub, in_images(i).im_lb, in_images(i).im_ub);                    
                end
                
            elseif isempty(option) || strcmp(option, 'single')
                for i=1:n
                    V = double(in_images(i).V); % convert to double precision
                    temp = obj.evaluate(V(:,:,:,1));
                    V(:,:,:,1) = temp;
                    images(i) = ImageStar(V, in_images(i).C, in_images(i).d, in_images(i).pred_lb, in_images(i).pred_ub, in_images(i).im_lb, in_images(i).im_ub);                    
                end
            else
                error('Unknown computation option');
            end
            
                      
        end
                 
    end
    
    
    methods(Static)
         % parse a trained input image layer from matlab
        function L = parse(input_image_layer)
            % @input_image_layer: input layer
            % @L: constructed layer
                        
            % author: Dung Tran
            % date: 8/1/2019
            
            
            if ~isa(input_image_layer, 'nnet.cnn.layer.ImageInputLayer')
                error('Input is not a Matlab nnet.cnn.layer.ImageInputLayer class');
            end
            
            L = ImageInputLayer(input_image_layer.Name, input_image_layer.InputSize, input_image_layer.AverageImage);
            fprintf('\nParsing a Matlab image input layer is done successfully');
            
        end
        
    end
    
    
    
end

