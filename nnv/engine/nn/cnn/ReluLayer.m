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
        
                
        
        % parse a trained relu Layer from matlab
        function parse(obj, relu_layer)
            % @relu_layer: a average pooling 2d layer from matlab deep
            % neural network tool box
                        
            % author: Dung Tran
            % date: 6/17/2019
            
            
            if ~isa(relu_Layer, 'reluLayer')
                error('Input is not a Matlab reluLayer class');
            end
            
            obj.Name = relu_layer.Name;
            fprintf('\nParsing a Matlab relu layer is done successfully');
            
        end
        
        
    end
        
    % exact reachability analysis using star set
    methods
        
        function images = reach_star_exact(~,in_image)
            % @in_image: an input imagestar
            % @image: output set
            
            % author: Dung Tran
            % date: 6/26/2019
             
            if ~isa(in_image, 'ImageStar')
                error('Input is not an imagestar');
            end
            
            
            nc = in_image.numChannel;
            h = in_image.height;
            w = in_image.width;
            np = in_image.numPred;
            
            N = h*w*nc; % total number of pixels in the input image         
            V(:, np+1) = zeros(N, 1);
            
            for i=1:np+1
                V(:,i) = reshape(in_image.V(:,:,i, :), N, 1);
            end
            
            I = Star(V, in_image.C, in_image.d, in_image.pred_lb, in_image.pred_ub);
            
            Y = PosLin.reach_star_exact(I,[]);
            n = length(Y);
            images(n) = ImageStar;
            for i=1:n
                images(i) = ImageStar(reshape(Y(i).V, [h, w, np+1, nc]), Y(i).C, Y(i).d, Y(i).predicate_lb, Y(i).predicate_ub);
            end
           
        end
        
        % reachability analysis for reluLayer using approximate imagestar
        function image = reach_star_approx(~,in_image)
            % @in_image: an input imagestar
            % @image: output set
            
            % author: Dung Tran
            % date: 6/27/2019
            
            if ~isa(in_image, 'ImageStar')
                error('Input is not an imagestar');
            end
            
            
            nc = in_image.numChannel;
            h = in_image.height;
            w = in_image.width;
            np = in_image.numPred;
            
            N = h*w*nc; % total number of pixels in the input image         
            V(:, np+1) = zeros(N, 1);
            
            for i=1:np+1
                V(:,i) = reshape(in_image.V(:,:,i, :), N, 1);
            end
            
            I = Star(V, in_image.C, in_image.d, in_image.pred_lb, in_image.pred_ub);
            Y = PosLin.reach_star_approx(I);
            image = ImageStar(reshape(Y.V, [h, w, Y.nVar+1, nc]), Y.C, Y.d, Y.predicate_lb, Y.predicate_ub);   

        end
        
            
    end
    
    
    
    
end

