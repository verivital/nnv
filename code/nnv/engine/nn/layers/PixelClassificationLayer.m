classdef PixelClassificationLayer < handle
    % Pixel Classification Layer object to verify segmentation networks
    % Author: Dung Tran
    % Date: 4/12/2020
    
    properties
        Name = 'PixelClassificationLayer';
        Classes = [];
        OutputSize = [];
    end
    
    methods
        
        function obj = PixelClassificationLayer(name, classes, outputSize)
            % @name: name of the layer
            % @classes: array of class
            % @classWeights: class weight used only for training process,
            % not for verification
            % @outputSize: outputSize
            % @lossFunction: lossfunction, we don't use it for verification
            
            % author: Dung Tran
            % date:4/12/2020
            
            if ~ischar(name)
                error('Invalid name, should be a charracter array');
            end
            
            if ~ismatrix(classes)
                error('Invalid classes, should be a matrix');
            end
                       
            if ~ismatrix(outputSize)
                error('Invalid outputSize');
            end           
            
            if length(outputSize) ~= 3
                error('Invalid outputSize matrix');
            end
            
            if length(classes) ~= outputSize(1, 3)
                error('Inconsistency betwen the number of classes and the outputSize');
            end
            
            obj.Name = name;
            obj.Classes = classes;
            obj.OutputSize = outputSize;
            
        end
            
    end
        
        
    methods
        
        % classified label for all pixels of an image
        function seg_im = evaluate(obj,image)
            %@image: an output image before softmax layer with the size of
            %        m1 x m2 x n, where n is the number of labels needs to
            %        be classified
            % @seg_im: segmentation image
            
            % author: Dung Tran
            % date: 4/12/2020
            
            n = size(image);
            if length(n)~= 3
                error('Output image should be a 3-dimensional image');
            end
            [~, y] = max(image, [], 3); 
            S = cell(n(1),n(2));
            classes = categories(obj.Classes);
           
            for i=1:n(1)
                for j=1:n(2)
                    S{i,j} =  classes{y(i,j)};
                end
            end
            
            seg_im = categorical(S);
        end
        
        % reachability with imagestar
        function seg_im = reach(obj, IS, method, reachOption)
            % @IS: imageStar input set
            % seg_im: segmentation image
            
            % author: Dung Tran
            % date: 4/12/2020
            
            [im_lb, im_ub] = IS.estimateRanges;
            [~,y1] = max(im_lb, [], 3);
            [~,y2] = max(im_ub, [], 3);
            classes = categories(obj.Classes);            
            n = size(y1);
            
            S = cell(n(1),n(2));
            for i=1:n(1)
                for j=1:n(2)
                    if y1(i,j) == y2(i,j)
                        S{i, j} = classes{y1(i,j)};
                    else
                        S{i,j} = 'unknown';
                    end
                end
            end
            
            seg_im = categorical(S);
        end
        
        
        
    end
    
    
    methods(Static)
        % parsing method
        
        function L = parse(pixel_classification_layer)
            % @pixel_classification_layer: 
            % @L: constructed layer
                        
            % author: Dung Tran
            % date: 4/12/2020
            
            
            if ~isa(pixel_classification_layer, 'nnet.cnn.layer.PixelClassificationLayer')
                error('Input is not a Matlab nnet.cnn.layer.PixelClassificationLayer');
            end
            
            L = PixelClassificationLayer(pixel_classification_layer.Name, pixel_classification_layer.Classes, pixel_classification_layer.OutputSize);
            
            fprintf('\nParsing a Matlab pixel classification layer is done successfully');
            
        end
        
        
    end
end

