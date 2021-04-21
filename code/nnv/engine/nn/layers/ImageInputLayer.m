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
        
    
    methods % reachability methods
        
        function image = reach_star_single_input(obj, in_image)
            % @in_image: an input ImageStar
            % @image: an output ImageStar
            
            % author: Dung Tran
            % date: 1/7/2020
            
            if ~isa(in_image, 'ImageStar')
                error('Input is not an ImageStar');
            end
            image = in_image.affineMap([], -obj.Mean);
        end
        
        % handling multiple inputs
        function images = reach_star_multipleInputs(obj, in_images, option)
            % @in_images: an array of ImageStars
            % @method: = 'exact-star' or 'approx-star' or 'abs-dom'
            % @option: = 'parallel' or 'single' or empty
            % @images: an array of ImageStar (if we use 'exact-star' method)
            %         or a single ImageStar set
            
            % author: Dung Tran
            % date: 1/7/2020
            
            n = length(in_images);
            images(n) = ImageStar;            
            if strcmp(option, 'parallel')
                parfor i=1:n
                    images(i) = obj.reach_star_single_input(in_images(i));
                end
            elseif strcmp(option, 'single') || isempty(option)
                for i=1:n
                    images(i) = obj.reach_star_single_input(in_images(i));
                end
            else
                error('Unknown computation option');

            end
        end
        
        function image = reach_zono(obj, in_image)
            % @in_image: an input ImageZono
            % @image: an output ImageZono
            
            % author: Dung Tran
            % date: 1/7/2020
            
            if ~isa(in_image, 'ImageZono')
                error('Input is not an ImageZono');
            end
            image = in_image.affineMap([], -obj.Mean);
        end
        
        % handling multiple inputs
        function images = reach_zono_multipleInputs(obj, in_images, option)
            % @in_images: an array of ImageZonos
            % @option: = 'parallel' or 'single' or empty
            % @images: an array of ImageZono 
            
            % author: Dung Tran
            % date: 1/7/2020
            
            n = length(in_images);
            images(n) = ImageZono;            
            if strcmp(option, 'parallel')
                parfor i=1:n
                    images(i) = obj.reach_zono(in_images(i));
                end
            elseif strcmp(option, 'single') || isempty(option)
                for i=1:n
                    images(i) = obj.reach_zono(in_images(i));
                end
            else
                error('Unknown computation option');

            end
        end
        
        
        
        function images = reach(varargin)
            % @in_image: an input imagestar
            % @image: output set
            % @option: = 'single' or 'parallel' 
            
            % author: Dung Tran
            % date: 6/26/2019
             
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
                    option = varargin{4};
                
                case 3
                    obj = varargin{1};
                    in_images = varargin{2};
                    method = varargin{3};
                    option = 'single';
                
                otherwise
                    error('Invalid number of input arguments (should be 2, 3, 4, 5, or 6)');
            end      
      
            if strcmp(method, 'approx-star') || strcmp(method, 'exact-star') || strcmp(method, 'abs-dom') || contains(method, "relax-star")
                images = obj.reach_star_multipleInputs(in_images, option);
            elseif strcmp(method, 'approx-zono')
                images = obj.reach_zono_multipleInputs(in_images, option);
            else
                error('Unknown reachability method');
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

