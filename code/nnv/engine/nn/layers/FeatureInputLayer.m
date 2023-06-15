classdef FeatureInputLayer < handle
    % The feature input layer class
    % Contain constructor and reachability analysis methods   
    % see: https://www.mathworks.com/help/deeplearning/ref/nnet.cnn.layer.featureinputlayer.html
    % 
    % Author: Diego Manzanas Lopez
    % Date: 03/15/2023
    
    properties % same properties as MATLAB's layer
        Name = 'FeatureInputLayer';
        InputSize = [];
        Normalization = 'none'; % default
        NormalizationDimension = 'auto'; % based on dimension of norm values
        Mean = []; 
        StandardDeviation = [];
        Min = [];
        Max = [];
    end
    
    
    methods % setting hyperparameters method
        
        % constructor of the class
        function obj = FeatureInputLayer(varargin)           
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
        
    
    methods % evaluation and reachbility methods
        
        % evaluate layer 
        function y = evaluate(obj, x)
            % @x: input image
            % @y: output image with normalization    

            if strcmp(obj.Normalization, 'none')
                y = double(x);
            elseif strcmp(obj.Normalization, 'zerocenter')
                y = double(x) - double(obj.Mean); % zerocenter nomalization
            else
                error('The normalization method is not supported yet.')
            end
                   
        end
        
        % Main star reachability function (TODO: update to other normalization techniques)
        function out_set = reach_star_single_input(obj, in_set)
            % @in_set: an input ImageStar
            % @out_set: an output ImageStar
            
            % Compute normalization
            if strcmp(obj.Normalization, 'none')
                out_set = in_set;
            elseif strcmp(obj.Normalization, 'zerocenter')
                out_set = in_set.affineMap([], -obj.Mean);
            else
                error('The normalization method is not supported yet.')
            end
            
        end
        
        % handling multiple inputs (star)
        function out_sets = reach_star_multipleInputs(obj, in_sets, option)
            % @in_sets: an array of ImageStars or Stars
            % @option: = 'parallel' or 'single' or empty
            % @out_sets: an array of ImageStar (if we use 'exact-star' method) or a single ImageStar set
            
            % check inputs
            setType = class(in_sets);
            n = length(in_sets);
            if ~contains(setType, 'Star')
                error('Wrong type of input sets to feature input layer. It must be a Star-based set.')
            elseif isa(in_sets, 'ImageStar')
                out_sets(n) = ImageStar;
            else
                out_sets(n) = Star;
            end
            
            % Begin computation
            if strcmp(option, 'parallel')
                parfor i=1:n
                    out_sets(i) = obj.reach_star_single_input(in_sets(i));
                end
            elseif strcmp(option, 'single') || isempty(option)
                for i=1:n
                    out_sets(i) = obj.reach_star_single_input(in_sets(i));
                end
            else
                error('Unknown computation option');

            end
        end
        
        % reachability with zonotopes
        function out_set = reach_zono(obj, in_set)
            % Operate on input set based on Normalization method
            if strcmp(obj.Normalization, 'none')
                out_set = in_set;
            elseif strcmp(obj.Normalization, 'zerocenter')
                out_set = in_set.affineMap([], -obj.Mean);
            else
                error('The normalization method is not supported yet.')
            end
        end
        
        % handling multiple zonotope inputs
        function images = reach_zono_multipleInputs(obj, in_images, option)
            
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
        
        % general reach function
        function images = reach(varargin)
             
            switch nargin
                
                 case 7
                    obj = varargin{1};
                    in_images = varargin{2};
                    method = varargin{3};
                    option = varargin{4};
                
                case 6
                    obj = varargin{1};
                    in_images = varargin{2};
                    method = varargin{3};
                    option = varargin{4};
                
                case 5
                    obj = varargin{1};
                    in_images = varargin{2};
                    method = varargin{3};
                    option = varargin{4};
                
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
    
    
    methods(Static) % parse method

         % parse a trained input feature layer from matlab
        function L = parse(layer)
            % @layer: input layer
            % @L: constructed layer

            if ~isa(layer, 'nnet.cnn.layer.FeatureInputLayer')
                error('Input is not a Matlab nnet.cnn.layer.FeatureInputLayer class');
            end
            L = FeatureInputLayer(layer.Name, layer.InputSize, layer.Normalization, layer.NormalizationDimension,...
                layer.Mean, layer.StandardDeviation, layer.Min, layer.Max);
        end
        
    end
    
end

