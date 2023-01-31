classdef SaturatingLinearLayer < ActivationFunctionLayer
    % The SaturationLinearLayer class in NN
    %   Contain constructor and reachability analysis methods
    %   Diego Manzanas: 12/06/2022
    %
    %   SATLIN : class for computing reachable set of Satlin Transfer Function 
    %      Reference: https://www.mathworks.com/help/deeplearning/ref/satlin.html
    % This looks equivalent to CLippedReluLayer, we could create this layer
    % parsing from there
    
    % setting hyperparameters method
    methods
        
        % constructor of the class
        function obj = SaturatingLinearLayer(varargin)           
            % author: Diego Manzanas
            % date: 12/06/2022   
            obj = obj@ActivationFunctionLayer(varargin);
        end
        
    end
        
    % evaluation method
    methods
        
        function y = evaluate(~, input)
            % @input: 2 or 3-dimensional array, for example, input(:, :, :), 
            % @y: 2 or 3-dimensional array, for example, y(:, :, :)
            
            n = size(input);
            N = 1;
            for i=1:length(n)
                N = N*n(i);
            end
            
            I = reshape(input, [N 1]);
            y = SatLin.evaluate(I);
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
            
            if ~isa(in_image, 'ImageStar') && ~isa(in_image, 'Star')
                error('input is not an ImageStar or Star');
            end
            
            if isa(in_image, "ImageStar")
                h = in_image.height;
                w = in_image.width;
                c = in_image.numChannel;
                            
                Y = SatLin.reach(in_image.toStar, method, [], relaxFactor, dis_opt, lp_solver); % reachable set computation with ReLU
                n = length(Y);
                images(n) = ImageStar;
                % transform back to ImageStar
                for i=1:n
                    images(i) = Y(i).toImageStar(h,w,c);
                end
            else
                images = SatLin.reach(in_image, method, [], relaxFactor, dis_opt, lp_solver); % reachable set computation with SatLin
            end

        end
        
        % reachability using ImageZono
        function image = reach_zono(~, in_image)
            % @in_image: an ImageZono input set
            
            if ~isa(in_image, 'ImageZono') && ~isa(in_image, "Zono")
                error('input is not an ImageZono or Zono');
            end
            
            if isa(in_image, "ImageZono")
                h = in_image.height;
                w = in_image.width;
                c = in_image.numChannels;
                In = in_image.toZono;
                Y = SatLin.reach(In, 'approx-zono');
                image = Y.toImageZono(h,w,c);
            else
                image = SatLin.reach(in_image, 'approx-zono');
            end
        end
                 
    end
    
end

