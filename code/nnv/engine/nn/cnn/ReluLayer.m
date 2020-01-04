classdef ReluLayer < handle
    % The Relu layer class in CNN
    %   Contain constructor and reachability analysis methods
    % Main references:
    % 1) An intuitive explanation of convolutional neural networks: 
    %    https://ujjwalkarn.me/2016/08/11/intuitive-explanation-convnets/
    % 2) More detail about mathematical background of CNN
    %    http://cs231n.github.io/convolutional-networks/
    %    http://cs231n.github.io/convolutional-networks/#pool
    % 3) Matlab implementation of Convolution2DLayer and MaxPooling (for training and evaluating purpose)
    %    https://www.mathworks.com/help/deeplearning/ug/layers-of-a-convolutional-neural-network.html
    %    https://www.mathworks.com/help/deeplearning/ref/nnet.cnn.layer.relulayer.html
    
    %   Dung Tran: 6/20/2019
    
    properties
        Name = 'relu_layer';        
    end
    
    
    % setting hyperparameters method
    methods
        
        % constructor of the class
        function obj = ReluLayer(varargin)           
            % author: Dung Tran
            % date: 6/26/2019    
            % update: 
            
            switch nargin
                
                case 1
                    
                    name = varargin{1};
                                        
                    if ~ischar(name)
                        error('Name is not char');
                    else
                        obj.Name = name;
                    end                    
                    
                case 0
                    
                    obj.Name = 'relu_layer';
                           
                otherwise
                    error('Invalid number of inputs (should be 0 or 1)');
                                 
            end 
             
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
            
            % author: Dung Tran
            % date: 6/26/2019
            
            
            n = size(input);
            N = 1;
            for i=1:length(n)
                N = N*n(i);
            end
            
            I = reshape(input, [N 1]);
            y = PosLin.evaluate(I);
            y = reshape(y, n);
                   
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
                    method = varargin{3};
                    option = [];
                    
                case 2
                    obj = varargin{1};
                    in_images = varargin{2};
                    method = [];
                    option = [];
                otherwise
                    error('Invalid number of input arguments (should be 1, 2 or 3)');
            end
      
            if isa(in_images(1), 'ImageStar') && ~isempty(method) && strcmp(method, 'approx-zono')
                error('Cannot use approx-zono reachability method for ImageStar input set, use exact-star, approx-star or abs-dom methods');
            end
            
            if isa(in_images(1), 'ImageZono') && ~isempty(method) && ~strcmp(method, 'approx-zono')
                error('Please use approx-zono method for the ImageZono input set, or transform ImageZono to ImageStar to perform the reachability analysis');
            end
            
            
            if isa(in_images(1), 'ImageStar') && strcmp(method, 'abs-dom')
                images = obj.reach_abs_dom(in_images, option); % use abstract-domain from DeepPoly
            else
                images = obj.reach_star(in_images, method, option); % use star method
            end
            
            if isa(in_images(1), 'ImageZono')
                images = obj.reach_zono(in_images, option); % use zonotope method from DeepZono
            end

        end
        
        
        % rechability using ImageStar
        function images = reach_star(~, in_images, method, option)
            % @in_images: an array of input ImageStar
            % @images: output images
            % @method: 'exact-star' or 'approx-star'
            % @option: 'parallel' or 'single': parallel computing
            
            % author: Dung Tran
            % date: 2/4/2020
            
            n = length(in_images);
            In(n) = Star; % preallocation
            h = in_images(1).height;
            w = in_images(1).width;
            c = in_images(1).numChannel;
            for i=1:n               
                if ~isa(in_images(i), 'ImageStar')
                    error('Input %d is not an imagestar', i);
                end
                In(i) = in_images(i).toStar;
            end
            
            Y = PosLin.reach(In, method, option); % reachable set computation with ReLU
            n = length(Y);
            images(n) = ImageStar;
            % transform back to ImageStar
            for i=1:n
                images(i) = Y(i).toImageStar(h,w,c);
            end
           
        end
        
        
        % reachability using ImageStar + Deepoly abstract-domain method
        function images = reach_abs_dom(~, in_images, option)
            % @in_images: an array of input ImageStar
            % @images: output images
            % @method: 'exact-star' or 'approx-star'
            % @option: 'parallel' or 'single': parallel computing
            
            % author: Dung Tran
            % date: 2/4/2020
            
            n = length(in_images);
            In(n) = Star; % preallocation
            h = in_images(1).height;
            w = in_images(1).width;
            c = in_images(1).numChannel;
            for i=1:n               
                if ~isa(in_images(i), 'ImageStar')
                    error('Input %d is not an imagestar', i);
                end
                In(i) = in_images(i).toStar;
            end
            
            images(n) = ImageStar;            
            if strcmp(option, 'parallel')
                parfor i=1:n
                    Y = PosLin.reach(In(i), 'abs-dom');
                    images(i) = Y.toImageStar(h,w,c);
                end
            elseif strcmp(option, 'single') || isempty(option)
                for i=1:n
                    Y = PosLin.reach(In(i), 'abs-dom');
                    images(i) = Y.toImageStar(h,w,c);
                end
            else
                error('Unknown computation option, should be parallel or single');
            end
           
        end
        
        
        % reachability using ImageZono
        function images = reach_zono(~,in_images, option)
            % @in_images: an array of input ImageZono
            % @images: output images
            % @option: 'parallel' or 'single' : parallel computing
            
            % author: Dung Tran
            % date: 2/4/2020
            
            n = length(in_images);
            In(n) = Zono; % preallocation
            h = in_images(1).height;
            w = in_images(1).width;
            c = in_images(1).numChannels;
            
            for i=1:n
                if ~isa(in_images(i), 'ImageZono')
                    error('Input %d is not an ImageZono', i);
                end                
                In(i) = in_images(i).toZono;                
            end
            
            images(n) = ImageZono;            
            if strcmp(option, 'parallel')
                parfor i=1:n
                    Y = PosLin.reach(In(i), 'approx-zono');
                    images(i) = Y.toImageZono(h,w,c);
                end
            elseif strcmp(option, 'single') || isempty(option)
                for i=1:n
                    Y = PosLin.reach(In(i), 'approx-zono');
                    images(i) = Y.toImageZono(h,w,c);
                end
            else
                error('Unknown computation option, should be parallel or single');
            end
            
        end
        
                 
    end
    
    
    methods(Static)
         % parse a trained relu Layer from matlab
        function L = parse(relu_layer)
            % @relu_layer: a average pooling 2d layer from matlab deep
            % neural network tool box
                        
            % author: Dung Tran
            % date: 6/17/2019
            
            
            if ~isa(relu_layer, 'nnet.cnn.layer.ReLULayer')
                error('Input is not a Matlab nnet.cnn.layer.ReLULayer class');
            end
            
            L = ReluLayer(relu_layer.Name);
            fprintf('\nParsing a Matlab relu layer is done successfully');
            
        end
        
    end
    
    
    
end

