classdef ElementwiseAffineLayer < handle
    % The ElementwiseAffineLayer layer class in CNN
    % author: Neelanjana Pal
    % date: 6/28/2021
    % update: Diego Manzanas
    %              10/31/2022
    %              Refactor properies and evaluation to mirror nnet.onnx.layer.ElementwiseAffineLayer
    
    properties
        Name = 'elementwise_affine_layer';
        Scale
        Offset
        DoScale
        DoOffset
    end
    
    
    % constructor
    methods
        
        % constructor of the class
        function obj = ElementwiseAffineLayer(varargin)           
            % author: Neelanjana Pal
            % date: 6/28/2021     
            
            switch nargin
                
                case 5
                    obj.Name = varargin{1};
                    obj.Scale = varargin{2};
                    obj.Offset = varargin{3};
                    obj.DoScale = varargin{4};
                    obj.DoOffset = varargin{5};
                case 4
                    obj.Scale = varargin{1};
                    obj.Offset = varargin{2};
                    obj.DoScale = varargin{3};
                    obj.DoOffset = varargin{4};
                otherwise
                    error('Invalid number of inputs (should be 4 or 5)');
            end 
             
        end
        
    end
        
    % evaluation method
    methods
        
        function y = evaluate(obj, x)
            % evaluate elementwise affine layer
            y = x;
            if obj.DoScale
                    y = y.*obj.Scale;
            end
            if obj.DoOffset
                y = y + obj.Offset; %reshape(obj.Offset, size(y));
            end          
        end 
       
    end   
     
 methods % reachability method
        
        %(reachability analysis using imagestar)
        function image = reach_star_single_input(obj, in_image)
            % @in_image: input imagestar
            % @image: output set
            
            if ~isa(in_image, 'ImageStar') && ~isa(in_image, "Star")
                error('Input set is not an ImageStar or Star');
            end
            
            if isa(in_image, "ImageStar")
                V = in_image.V; % Initialize value dimensions
    %             n = in_image.numPred;
                % Process scaling first
                if obj.DoScale
                    if isscalar(obj.Scale)
                        V = double(obj.Scale)*in_image.V;
                    elseif length(size(obj.Scale)) == length(size(in_image.V)) && size(obj.Scale) == size(in_image.V)
                        V = double(obj.Scale).*in_image;
                    elseif length(obj.Scale) == size(in_image.V,1)
                        for k=1:length(obj.Scale)
                            V(k, : , : , :) = V(k, : , : , :) * obj.Scale(k);
                        end
                    elseif length(obj.Scale) == size(in_image.V,2)
                        for k=1:length(obj.Scale)
                            V(:, k , : , :) = V(:, k , : , :) * obj.Scale(k);
                        end
                    elseif length(obj.Scale) == size(in_image.V,3)
                        for k=1:length(obj.Scale)
                            V(:, : , k , :) = V(:, : , k , :) * obj.Scale(k);
                        end
                    elseif length(obj.Scale) == size(in_image.V,4)
                        for k=1:length(obj.Scale)
                            V(:, : , : , k) = V(:, : , : , k) * obj.Scale(k);
                        end
                    else
                        error('TODO: add support for other tensor shapes (Scale)')
                    end
                end
    
                % Process offset last
                if obj.DoOffset
                    offset = squeeze(obj.Offset);
    %                 a = size(offset);
                    if isscalar(obj.Offset)
                        V = V + obj.Offset;
                    elseif length(obj.Offset) == size(in_image.V,1)
                        for k=1:length(obj.Offset)
                            V(k, : , : , :) = V(k, : , : , :) + obj.Offset(k);
                        end
                    elseif length(obj.Offset) == size(in_image.V,2)
                        for k=1:length(obj.Offset)
                            V(:, k , : , :) = V(:, k , : , :) + obj.Offset(k);
                        end
                    elseif length(obj.Offset) == size(in_image.V,3)
                        for k=1:length(obj.Offset)
                            V(:, : , k , 1) = V(:, : , k , 1) + obj.Offset(k); % bias only added to first column
                        end
                    elseif length(obj.Offset) == size(in_image.V,4)
                        for k=1:length(obj.Offset)
                            V(:, : , : , k) = V(:, : , : , k) + obj.Offset(k);
                        end
                    else
                        error('TODO: add support for other tensor shapes (Offset).')
                    end
                end
    
                % return output set
                image = ImageStar(V, in_image.C, in_image.d, in_image.pred_lb, in_image.pred_ub);
            else % reachability with Star sets
                % In general, this layer is used as the bias addition following a FullyConnectedLayer 
                % Scale (weights)
                image = in_image; % create copy to perform operations (if need to)
                if obj.DoScale
                    image = image.affineMap(obj.Scale, []); % W*x
                end
                % Offset (bias)
                if obj.DoOffset
                    image = image.affineMap(diag(ones(1,a.dim)), obj.Offset); % x + b
                end
            end
            
        end
        
        % handle multiple inputs
        function S = reach_star_multipleInputs(obj, inputs, option)
            % @inputs: an array of ImageStars
            % @option: = 'parallel' or 'single'
            % @S: output ImageStar
            
            n = length(inputs);
            S(n) = ImageStar;
            if strcmp(option, 'parallel')
                parfor i=1:n
                    S(i) = obj.reach_star_single_input(inputs(i));
                end
            elseif strcmp(option, 'single') || isempty(option)
                for i=1:n
                    S(i) = obj.reach_star_single_input(inputs(i));
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
            
            if obj.DoScale || obj.DoOffset
                if strcmp(method, 'approx-star') || strcmp(method, 'exact-star') || contains(method, "relax-star")
                    IS = obj.reach_star_multipleInputs(in_images, option);
                else
                    error('Unsupported/Unknown reachability method, only star based methods are supported');
                end
            else
                IS = in_images;
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
                        
            L = ElementwiseAffineLayer(elementwise_affine_layer.Name, elementwise_affine_layer.Scale, elementwise_affine_layer.Offset, ...
                elementwise_affine_layer.DoScale, elementwise_affine_layer.DoOffset);         
            
        end
        
        
    end
    
    
    
    
end

