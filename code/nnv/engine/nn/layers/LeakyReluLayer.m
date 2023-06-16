classdef LeakyReluLayer < ActivationFunctionLayer
    % The LeakyRelu layer class in NN
    %   Contain constructor and reachability analysis methods
    
    %   Diego Manzanas: 12/06/2022
    %
    
    properties
        gamma = 0.01; % MATLAB default scale value
    end
    
    % setting hyperparameters method (custom to account for gamma)
    methods
        
        % constructor of the class
        function obj = LeakyReluLayer(varargin)           
            % author: Diego Manzanas
            % date: 12/06/2022
            
            obj = obj@ActivationFunctionLayer(varargin(1:5));
            obj.gamma = varargin{6}; % scale in MATLAB LeakyRelu layer
             
        end
        
    end
        
    % evaluation method
    methods
        
        function y = evaluate(obj, input)
            % @input: 2 or 3-dimensional array, for example, input(:, :, :), 
            % @y: 2 or 3-dimensional array, for example, y(:, :, :)
            
            n = size(input);
            N = 1;
            for i=1:length(n)
                N = N*n(i);
            end
            
            I = reshape(input, [N 1]);
            y = LeakyReLU.evaluate(I, obj.gamma);
            y = reshape(y, n);
                   
        end
        
		function y = evaluateSequence(obj, input)
            % @input: 2 or 3-dimensional array, for example, input(:, :, :), 
            % @y: 2 or 3-dimensional array, for example, y(:, :, :)
            
            % author: Neelanjana Pal
            % date: 1/6/2023
            
            y = LeakyReLU.evaluate(input, obj.gamma);
                   
        end										 
    end
    
    methods % reachability methods
        
        % reachability using ImageStar
        function images = reach_star_single_input(obj, in_image, method, relaxFactor, dis_opt, lp_solver)
            % @in_image: an ImageStar input set
            % @method: = 'exact-star' or 'approx-star' or 'abs-dom'
            % @relaxFactor: of approx-star method
            % @dis_opt: display option = [] or 'display'
            % @lp_solver: lp solver
            % @images: an array of ImageStar (if we use 'exact-star' method)
            %         or a single ImageStar set
                        
            if ~isa(in_image, 'ImageStar') && ~isa(in_image, 'Star')
                error('input is not an ImageStar or Star');
            end
            option = [];
            if isa(in_image, "ImageStar")
                h = in_image.height;
                w = in_image.width;
                c = in_image.numChannel;
                
                Y = LeakyReLU.reach(in_image.toStar, obj.gamma, method, option, relaxFactor, dis_opt, lp_solver); % reachable set computation with LeakyReLU
                n = length(Y);
                images(n) = ImageStar;
                % transform back to ImageStar
                for i=1:n
                    images(i) = Y(i).toImageStar(h,w,c);
                end
            else
                images = LeakyReLU.reach(in_image, obj.gamma, method, option, relaxFactor, dis_opt, lp_solver); % reachable set computation with LeakyReLU
            end
        end
        
        % reachability using ImageZono
        function image = reach_zono(obj, in_image)
            % @in_image: an ImageZono input set
            
            if ~isa(in_image, 'ImageZono') && ~isa(in_image, "Zono")
                error('input is not an ImageZono or Zono');
            end
            
            if isa(in_image, "ImageZono")
                h = in_image.height;
                w = in_image.width;
                c = in_image.numChannels;
                In = in_image.toZono;
                Y = LeakyReLU.reach(In, obj.gamma, 'approx-zono');
                image = Y.toImageZono(h,w,c);
            else
                image = LeakyReLU.reach(in_image, obj.gamma, 'approx-zono');
            end
            
        end
                 
    end
    
    
    methods(Static)
         % parse a trained leaky relu Layer from matlab
        function L = parse(layer)
            if ~isa(layer, 'nnet.cnn.layer.LeakyReLULayer')
                error('Input is not a Matlab nnet.cnn.layer.LeakyReLULayer class');
            end
            L = LeakyReluLayer(layer.Name, layer.NumInputs, layer.InputNames, layer.NumOutputs, layer.OutputNames, layer.Scale);
        end
        
    end
    
end

