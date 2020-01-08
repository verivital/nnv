classdef ImageInputLayer < handle
    % The Image input layer class in CNN
    %   Contain constructor and reachability analysis methods   
    %   Dung Tran: 8/1/2019
    %   update: 1/2/2020 (Dung Tran)
    %       update reason: different Matlab versions have different names of AverageImage
    %       for example: In 2018b: Matlab uses the name "AverageImage"
    %                    In 2019b: Matlab uses the name "Mean" 
    %       We decided to update the name AverageImage to Mean in 1/2/2020
    %       
    
    properties
        Name = 'ImageInputLayer';
        InputSize = [];
        Normalization = 'zerocenter'; %default
        Mean = []; % in 
                
    end
    
    
    % setting hyperparameters method
    methods
        
        % constructor of the class
        function obj = ImageInputLayer(varargin)           
            % author: Dung Tran
            % date: 6/26/2019    
            % update: 1/2/2020
            %   update reason: change the way the object receives input
            %   arguments
            
            
            if mod(nargin, 2) ~= 0
                error('Invalid number of input arguments');
            end
            
            for i=1:nargin-1
                
                if mod(i, 2) ~= 0
                    
                    if strcmp(varargin{i}, 'Name')
                        obj.Name = varargin{i+1};
                    elseif strcmp(varargin{i}, 'InputSize')
                        obj.InputSize = varargin{i+1};
                    elseif strcmp(varargin{i}, 'Mean')
                        obj.Mean = double(varargin{i+1});
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
            % date: 8/1/2019
                             
            if isempty(obj.Mean)
                y = double(input);
            else
                y = double(input) - double(obj.Mean); % zerocenter nomalization
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
        
        
        % reachability analysis using ImageStar or ImageZono method
        function images = reach_zono(varargin)
            % @in_images: an array of input ImageZono
            % @option: 'parallel' or 'single' or '[]'
            % @images: output set
            
            % author: Dung Tran
            % date: 1/4/2020
            
            
            switch nargin
                case 3
                    obj = varargin{1};
                    in_images = varargin{2};
                    option = varargin{3};
                case 2
                    obj = varargin{1};
                    in_images = varargin{2};
                    option = 'single';
                otherwise
                    error('Invalid number of inputs, should be 1 or 2');
            end
        
            
            
            n = length(in_images);
            for i=1:n
                if ~isa(in_images(i), 'ImageZono') || ~isa(in_images(i), 'ImageStar')
                    error('The %d the input is not an ImageStar or ImageZono');
                end                
            end
            
            if isa(in_images(1), 'ImageStar')
                images(n) = ImageStar;
            elseif isa(in_images(1), 'ImageZono')
                images(n) = ImageZono;
            end
   
            mean_image = obj.Mean;
            if strcmp(option, 'parallel')
                parfor i=1:n
                    images(i) = in_images(i).affineMap([], mean_image);
                end
            elseif strcmp(option, 'single') || isempty(option)
                for i=1:n
                    images(i) = in_images(i).affineMap([], mean_image);
                end
            else
                error('Unknown computation option, should be parallel or single');
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
            
            if isprop(input_image_layer, 'Mean')
                L = ImageInputLayer('Name', input_image_layer.Name, 'InputSize', input_image_layer.InputSize, 'Mean', input_image_layer.Mean);
            elseif isprop(input_image_layer, 'AverageImage')
                L = ImageInputLayer('Name', input_image_layer.Name, 'InputSize', input_image_layer.InputSize, 'Mean', input_image_layer.AverageImage);
            else
                error('Mean or AverageImage property does not exist in the Input Image Layer');
            end
            
            fprintf('\nParsing a Matlab image input layer is done successfully');
            
        end
        
    end
    
    
    
end

