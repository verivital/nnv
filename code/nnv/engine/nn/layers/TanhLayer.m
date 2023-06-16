classdef TanhLayer < ActivationFunctionLayer
    % The TanhdLayer class in NN
    %   Contain constructor and reachability analysis methods
    %   Dung Tran: 6/9/2020
    %   
    %   update: Diego Manzanas, 12/07/2022
    %          - Inheritance using ActivationFunctionLayer
    
    % setting hyperparameters method
    methods
        
        % constructor of the class
        function obj = TanhLayer(varargin)           
            % author: Dung Tran
            % date: 6/9/2020    
            %   
            % update: Diego Manzanas, 12/07/2022
            %   - Inheritance using ActivationFunctionLayer
            obj = obj@ActivationFunctionLayer(varargin);
        end
        
        
    end
        
    % evaluation method
    methods
        
        function y = evaluate(~, input)
            % @input: 2 or 3-dimensional array, for example, input(:, :, :), 
            % @y: 2 or 3-dimensional array, for example, y(:, :, :)
            
            % author: Dung Tran
            % date: 6/26/2019
            
            % @y: high-dimensional array (output volume)
            
            n = size(input);
            N = 1;
            for i=1:length(n)
                N = N*n(i);
            end
            
            I = reshape(input, [N 1]);
            y = TanSig.evaluate(I);
            y = reshape(y, n);
                   
        end
		function y = evaluateSequence(~, input)
            % @input: 2 or 3-dimensional array, for example, input(:, :, :), 
            % @y: 2 or 3-dimensional array, for example, y(:, :, :)
            
            % author: Neelanjana Pal
            % date: 1/6/2023
            
            y = TanSig.evaluate(input);
                   
        end									   
        
    end
        
    
    methods % reachability methods
        
        % reachability using ImageStar
        function images = reach_star_single_input(~, in_image, method, relaxFactor, dis_opt, lp_solver)
            % @in_image: an ImageStar input set
            % @method: = 'exact-star' or 'approx-star' or 'abs-dom'
            % @relaxFactor: for approx-star method only
            % @images: an array of ImageStar (if we use 'exact-star' method)
            %         or a single ImageStar set
            
            % author: Dung Tran
            % date: 6/9/2020
            % update: 6/26/2020: add relaxed approx-star method
            %         7/16/2020: add display option + lp_solver option
            
            if ~isa(in_image, 'ImageStar') && ~isa(in_image, 'Star')
                error('input is not an ImageStar or Star');
            end
            
            if isa(in_image, "ImageStar")
                h = in_image.height;
                w = in_image.width;
                c = in_image.numChannel;
                
                Y = TanSig.reach(in_image.toStar, method, [], relaxFactor, dis_opt, lp_solver); % reachable set computation with Tanh
                n = length(Y);
                images(n) = ImageStar;
                % transform back to ImageStar
                for i=1:n
                    images(i) = Y(i).toImageStar(h,w,c);
                end
            else
                images = TanSig.reach(in_image, method, [], relaxFactor, dis_opt, lp_solver); % reachable set computation with Tanh
            end

        end
        
        % reachability using ImageZono
        function image = reach_zono(~, in_image)
            % @in_image: an ImageZono input set
            
            % author: Dung Tran
            % date: 6/9/2020
            
            if ~isa(in_image, 'ImageZono') && ~isa(in_image, "Zono")
                error('input is not an ImageZono or Zono');
            end
            
            if isa(in_image, "ImageZono")
                h = in_image.height;
                w = in_image.width;
                c = in_image.numChannels;
                In = in_image.toZono;
                Y = TanSig.reach(In, 'approx-zono');
                image = Y.toImageZono(h,w,c);
            else
                image = TanSig.reach(in_image, 'approx-zono');
            end

        end
                 								   
    end
    
    
    methods(Static)
         % parse a trained tanh layer from matlab
        function L = parse(layer)
            % @layer: a tanh layer from Matlab or Keras (parsed from importKerasNetwork)
                        
            % author: Dung Tran
            % date: 6/9/2020          
            
            if ~isa(layer, 'nnet.keras.layer.TanhLayer') && ~isa(layer, 'nnet.cnn.layer.TanhLayer')  
                error('Input is not a sigmoid layer');
            end
            
            L = TanhLayer(layer.Name, layer.NumInputs, layer.InputNames, layer.NumOutputs, layer.OutputNames);
        end
        
    end
    
end

