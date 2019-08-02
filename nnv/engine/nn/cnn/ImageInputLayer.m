classdef ImageInputLayer < handle
    % The Image input layer class in CNN
    %   Contain constructor and reachability analysis methods   
    %   Dung Tran: 8/1/2019
    
    properties
        Name = 'ImageInputLayer';
        InputSize = [];
        Normalization = 'zerocenter'; %default
        AverageImage = []; % average image
        
    end
    
    
    % setting hyperparameters method
    methods
        
        % constructor of the class
        function obj = ImageInputLayer(varargin)           
            % author: Dung Tran
            % date: 6/26/2019    
            % update: 
            
            switch nargin
                
                case 2
                    
                    name = varargin{1};
                                                            
                    if ~ischar(name)
                        error('Name is not char');
                    else
                        obj.Name = name;
                    end                    
                    
                    inputSize = varargin{2};
                    obj.InputSize = inputSize;
                                        
                case 3
                    
                    name = varargin{1};
                    inputSize = varargin{2};
                    averageImage = varargin{3};
                                                            
                    if ~ischar(name)
                        error('Name is not char');
                    else
                        obj.Name = name;
                    end                    
                    
                    
                    obj.InputSize = inputSize;
                    
                    averageImage = reshape(averageImage, inputSize);
                    obj.AverageImage = averageImage;
               
                case 0
                    
                    obj.Name = 'InputImageLayer';
                           
                otherwise
                    error('Invalid number of inputs (should be 0, 2, or 3)');
                                 
            end 
             
        end
        
        
    end
        
    % evaluation method
    methods
        
        function y = evaluate(obj, input)
            % @input: input image
            % @y: output image with normalization
            
            % author: Dung Tran
            % date: 8/1/2019
                             
            if isempty(obj.AverageImage)
                y = double(input);
            else
                y = double(input) - double(obj.AverageImage); % zerocenter nomalization
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

