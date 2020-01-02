classdef BatchNormalizationLayer < handle
    % The Batch Normalization Layer class in CNN
    %   Contain constructor and reachability analysis methods   
    %   Dung Tran: 1/1/2020
    
    properties
        Name = 'BatchNormalizationLayer';
        
        NumChannels = [];
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
                
                if mod(i, 2) ~= 0
                    
                    if strcmp(varargin{i}, 'Name')
                        obj.Name = varargin{i+1};
                    elseif strcmp(varargin{i}, 'NumChannels')
                        obj.NumChannels = varargin{i+1};
                    elseif strcmp(varargin{i}, 'TrainedMean')
                        obj.TrainedMean = double(varargin{i+1});
                    elseif strcmp(varargin{i}, 'TrainedVariance')
                        obj.TrainedVariance = double(varargin{i+1});
                    elseif strcmp(varargin{i}, 'Epsilon')
                        obj.Epsilon = double(varargin{i+1});
                    elseif strcmp(varargin{i}, 'Offset')
                        obj.Offset = double(varargin{i+1});
                    elseif strcmp(varargin{i}, 'Scale')
                        obj.Scale = double(varargin{i+1});
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
                y = double(input);
                for i=1:obj.NumChannels
                    y(:,:,i) = y(:,:,i) - obj.TrainedMean(1,1,i);
                    y(:,:,i) = y(:,:,i)/(sqrt(obj.TrainedVariance(1,1,i) + obj.Epsilon));
                    y(:,:,i) = obj.Scale(1, 1, i)*y(:,:,i) + obj.Offset(1,1,i);
                end
                
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
            
            if isempty(obj.TrainedMean) || isempty(obj.TrainedVariance) || isempty(obj.Epsilon) || isempty(obj.Offset) || isempty(obj.Scale)
                error('Batch Normalization Layer does not have enough parameters');
            end
            
            if strcmp(option, 'parallel')
                parfor i=1:n
                    x = in_images(i).affineMap(1/sqrt(obj.TrainedVariance^2 + obj.Epsilon), []);
                    x = x.affineMap([], obj.TrainedMean);
                    images(i) = x.affineMap(obj.Scale, obj.Offset);
                end
                
            elseif isempty(option) || strcmp(option, 'single')
                for i=1:n
                    x = in_images(i).affineMap(1/sqrt(obj.TrainedVariance^2 + obj.Epsilon), []);
                    x = x.affineMap([], obj.TrainedMean);
                    images(i) = x.affineMap(obj.Scale, obj.Offset);
                end
            else
                error('Unknown computation option');
            end
            
                      
        end
                 
    end
    
    
    methods(Static)
         % parse a trained batch normalization layer from matlab
        function L = parse(layer)
            % @layer: batch normalization layer
            % @L: constructed layer
                        
            % author: Dung Tran
            % date: 1/1/2020
            
            
            if ~isa(layer, 'nnet.cnn.layer.BatchNormalizationLayer')
                error('Input is not a Matlab nnet.cnn.layer.BatchNormalizationLayer class');
            end
            
            display(layer.Name);
            display(layer.TrainedMean);
            display(layer.TrainedVariance);
            display(layer.NumChannels);
            
            L = BatchNormalizationLayer('Name', layer.Name, 'NumChannels', layer.NumChannels, 'TrainedMean', layer.TrainedMean, 'TrainedVariance', layer.TrainedVariance, 'Epsilon', layer.Epsilon, 'Offset', layer.Offset, 'Scale', layer.Scale);
            fprintf('\nParsing a Matlab batch normalization layer is done successfully');
            
        end
        
    end
    
    
    
end

