classdef SigmoidLayer < ActivationFunctionLayer
    % The SigmoidLayer class in NN
    %   Contain constructor and reachability analysis methods
    %   Dung Tran: 6/9/2020
    %   
    %   update: Diego Manzanas, 12/07/2022
    %          - Inheritance using ActivationFunctionLayer
    
    % setting hyperparameters method
    methods
        
        % constructor of the class
        function obj = SigmoidLayer(varargin)           
            % author: Dung Tran
            % date: 6/9/2020   
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
            y = LogSig.evaluate(I);
            y = reshape(y, n);
                   
        end
        
    end
        
    
    methods % reachability methods
        
        % reachability using ImageStar
        function images = reach_star_single_input(~, in_image, method, relaxFactor, dis_opt, lp_solver)
            % @in_image: an ImageStar input set
            % @method: = 'exact-star' or 'approx-star' or 'abs-dom'
            % @relaxFactor: for approx-star method
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
                            
                Y = LogSig.reach(in_image.toStar, method, [], relaxFactor, dis_opt, lp_solver); % reachable set computation with sigmoid (logsig)
                n = length(Y);
                images(n) = ImageStar;
                % transform back to ImageStar
                for i=1:n
                    images(i) = Y(i).toImageStar(h,w,c);
                end
            else
                images = LogSig.reach(in_image, method, [], relaxFactor, dis_opt, lp_solver); % reachable set computation with sigmoid (logsig)
            end

        end
        
        % reachability using ImageZono
        function image = reach_zono(~, in_image)
            % @in_image: an ImageZono input set
            
            % author: Dung Tran
            % date: 6/9/2020
            
            if ~isa(in_image, 'ImageZono') && ~isa(in_image, 'Zono')
                error('Input is not an ImageZono or Zono');
            end
            
            if isa(in_image, 'ImageZono')
                h = in_image.height;
                w = in_image.width;
                c = in_image.numChannels;
                In = in_image.toZono;
                Y = LogSig.reach(In, 'approx-zono');
                image = Y.toImageZono(h,w,c);
            else
                image = LogSig.reach(in_image, 'approx-zono');
            end
        end
                 
    end
    
    
    methods(Static)
         % parse a trained simoid layer from matlab
        function L = parse(layer)
            % @sigmoid_layer: 
                        
            % author: Dung Tran
            % date: 6/9/2020
            % modified: added onnx flatten sigmoid layers
            %       by: Neelanjana Pal
            % date: 6/25/2021
            if ~isa(layer, 'nnet.keras.layer.SigmoidLayer') && ~isa(layer,'nnet.onnx.layer.SigmoidLayer') && ~isa(layer,'nnet.cnn.layer.SigmoidLayer')
                error('Input is not a sigmoid layer');
            end
            
            L = SigmoidLayer(layer.Name, layer.NumInputs, layer.InputNames, layer.NumOutputs, layer.OutputNames);
            
        end
        
    end
    
end

