classdef Image3DInputLayer < handle
    % The Image 3d input layer 
    %   Contain constructor and reachability analysis methods   
    %   Diego Manzanas: 10/11/2023
    
    properties % same properties as MATLAB's layer
        Name = 'Image3dInputLayer';
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
        function obj = Image3DInputLayer(varargin)
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
            
            if ~isa(in_image, 'VolumeStar')
                error('Input is not an VolumeStar');
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
            
            n = length(in_images);
            images(n) = VolumeStar;

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
            
            if ~isa(layer, 'nnet.cnn.layer.Image3DInputLayer')
                error('Input is not a Matlab nnet.cnn.layer.Image3DInputLayer class');
            end
            
            L = Image3DInputLayer(layer.Name, layer.InputSize, layer.Normalization, layer.NormalizationDimension,...
                layer.Mean, layer.StandardDeviation, layer.Min, layer.Max);
                        
        end
        
    end
    
    
end

