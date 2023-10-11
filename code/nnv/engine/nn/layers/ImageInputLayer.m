classdef ImageInputLayer < handle
    % The Image input layer 
    %   Contain constructor and reachability analysis methods   
    %   Dung Tran: 8/1/2019
    %   update: 1/2/2020 (Dung Tran)
    %       update reason: different Matlab versions have different names of AverageImage
    %       for example: In 2018b: Matlab uses the name "AverageImage"
    %                    In 2019b: Matlab uses the name "Mean" 
    %       We decided to update the name AverageImage to Mean in 1/2/2020
    %
    %  update: 03/16/2023 (Diego Manzanas)
    %       Add properties and update routines for different normalization methods
    
    properties % same properties as MATLAB's layer
        Name = 'ImageInputLayer';
        InputSize = [];
        Normalization = 'none'; % default
        NormalizationDimension = 'auto'; % based on dimension of norm values
        Mean = []; 
        StandardDeviation = [];
        Min = [];
        Max = [];

    end
    
    
    % create layer
    methods
        
        % constructor of the class 
        function obj = ImageInputLayer(varargin)
            % assign values to properties
            switch nargin
                case 0 % create empty layer
                case 1 % define input size
                    obj.InputSize = varargin{1};
                case 2 % input size + normalization
                    obj.InputSize = varargin{1};
                    obj.Normalization = varargin{2};
                case 8 % define all properties
                    obj.Name = varargin{1};
                    obj.InputSize = varargin{2};
                    obj.Normalization = varargin{3};
                    obj.NormalizationDimension = varargin{4};
                    obj.Mean = varargin{5};
                    obj.StandardDeviation = varargin{6};
                    obj.Min = varargin{7};
                    obj.Max = varargin{8};
                otherwise
                    error("Wrong number of inputs, it must be 0,1,2, or 8.")
            end
        end
        
    end
        
    % evaluation method
    methods
        
        % TODO: add normalization options
        function y = evaluate(obj, x)
            % @x: input image
            % @y: output image with normalization
            
            if strcmp(obj.Normalization, 'none')
                y = x;
            elseif strcmp(obj.Normalization, 'zerocenter')
                y = x - obj.Mean; % zerocenter nomalization
            else
                error('The normalization method is not supported yet.')
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

            % Compute normalization
            if strcmp(obj.Normalization, 'none')
                image = in_image;
            elseif strcmp(obj.Normalization, 'zerocenter')
                image = in_image.affineMap([], -obj.Mean);
            else
                error('The normalization method is not supported yet.')
            end

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

            % Compute normalization
            if strcmp(obj.Normalization, 'none')
                image = in_image;
            elseif strcmp(obj.Normalization, 'zerocenter')
                image = in_image.affineMap([], -obj.Mean);
            else
                error('The normalization method is not supported yet.')
            end
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
        function L = parse(layer)
            % @layer: input layer
            % @L: constructed layer            
            
            if ~isa(layer, 'nnet.cnn.layer.ImageInputLayer')
                error('Input is not a Matlab nnet.cnn.layer.ImageInputLayer class');
            end
            
            L = ImageInputLayer(layer.Name, layer.InputSize, layer.Normalization, layer.NormalizationDimension,...
                layer.Mean, layer.StandardDeviation, layer.Min, layer.Max);
                        
        end
        
    end
    
    
    
end

