classdef TanhLayer < handle
    % The TanhdLayer class in CNN
    %   Contain constructor and reachability analysis methods
    %   Dung Tran: 6/9/2020
    
    properties
        Name = 'tanh_layer';
        
        NumInputs = 1;
        InputNames = {'in'};
        NumOutputs = 1;
        OutputNames = {'out'};
    end
    
    
    % setting hyperparameters method
    methods
        
        % constructor of the class
        function obj = TanhLayer(varargin)           
            % author: Dung Tran
            % date: 6/9/2020    
            % update: 
            
            switch nargin
                
                case 5 % used for parsin a matlab relu layer
                    obj.Name = varargin{1};
                    obj.NumInputs = varargin{2};
                    obj.InputNames = varargin{3};
                    obj.NumOutputs = varargin{4};
                    obj.OutputNames = varargin{5};

                case 1
                    
                    name = varargin{1};
                                        
                    if ~ischar(name)
                        error('Name is not char');
                    else
                        obj.Name = name;
                    end                    
                    
                case 0
                    
                    obj.Name = 'tanh_layer';
                           
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
            y = TanSig.evaluate(I);
            y = reshape(y, n);
                   
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
            
            if ~isa(in_image, 'ImageStar')
                error('input is not an ImageStar');
            end
            
            if relaxFactor < 0
                error('Invalid relaxFactor');
            end
            
            h = in_image.height;
            w = in_image.width;
            c = in_image.numChannel;
            
            Y = TanSig.reach(in_image.toStar, method, [], relaxFactor, dis_opt, lp_solver); % reachable set computation with ReLU
            n = length(Y);
            images(n) = ImageStar;
            % transform back to ImageStar
            for i=1:n
                images(i) = Y(i).toImageStar(h,w,c);
            end

        end
        
        
        % hangling multiple inputs
        function images = reach_star_multipleInputs(obj, in_images, method, option, relaxFactor, dis_opt, lp_solver)
            % @in_images: an array of ImageStars
            % @method: = 'exact-star' or 'approx-star' or 'abs-dom'
            % @option: = 'parallel' or 'single' or empty
            % @relaxFactor: for approx-star method
            % @images: an array of ImageStar (if we use 'exact-star' method)
            %         or a single ImageStar set
            
            % author: Dung Tran
            % date: 6/9/2020
            % update: 6/26/2020: add relaxed approx-star method
            %         7/16/2020: add display option + lp_solver option
            
            
            images = [];
            n = length(in_images);
                        
            if strcmp(option, 'parallel')
                parfor i=1:n
                    images = [images obj.reach_star_single_input(in_images(i), method, relaxFactor, dis_opt, lp_solver)];
                end
            elseif strcmp(option, 'single') || isempty(option)
                for i=1:n
                    images = [images obj.reach_star_single_input(in_images(i), method, relaxFactor, dis_opt, lp_solver)];
                end
            else
                error('Unknown computation option');

            end
            
        end
        
        
        % reachability using ImageZono
        function image = reach_zono(~, in_image)
            % @in_image: an ImageZono input set
            
            % author: Dung Tran
            % date: 6/9/2020
            
            if ~isa(in_image, 'ImageZono')
                error('input is not an ImageZono');
            end
            
            h = in_image.height;
            w = in_image.width;
            c = in_image.numChannels;
            In = in_image.toZono;
            Y = TanSig.reach(In, 'approx-zono');
            image = Y.toImageZono(h,w,c);
            
        end
        
        % handling multiple inputs
        function images = reach_zono_multipleInputs(obj, in_images, option)
            % @in_images: an array of ImageStars
            % @option: = 'parallel' or 'single' or empty
            % @images: an array of ImageStar (if we use 'exact-star' method)
            %         or a single ImageStar set
            
            % author: Dung Tran
            % date: 6/9/2020
            
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
        
        
        % MAIN REACHABILITY METHOD
        function images = reach(varargin)
            % @in_image: an input imagestar
            % @image: output set
            % @option: = 'single' or 'parallel' 
            
            % author: Dung Tran
            % date: 6/9/2020
            % update: 7/16/2020: add display option + lp_solver option
             
            switch nargin
                
                case 7
                    obj = varargin{1};
                    in_images = varargin{2};
                    method = varargin{3};
                    option = varargin{4};
                    relaxFactor = varargin{5}; % for relaxed approx-star method
                    dis_opt = varargin{6};
                    lp_solver = varargin{7};
                    
                case 6
                    obj = varargin{1};
                    in_images = varargin{2};
                    method = varargin{3};
                    option = varargin{4};
                    relaxFactor = varargin{5}; % for relaxed approx-star method
                    dis_opt = varargin{6};
                    lp_solver = 'glpk';
                case 5
                    obj = varargin{1};
                    in_images = varargin{2};
                    method = varargin{3};
                    option = varargin{4};
                    relaxFactor = varargin{5}; % for relaxed approx-star method
                    dis_opt = [];
                    lp_solver = 'glpk';
                case 4
                    obj = varargin{1};
                    in_images = varargin{2};
                    method = varargin{3};
                    option = varargin{4};
                    relaxFactor = 0;
                    dis_opt = [];
                    lp_solver = 'glpk';
                
                case 3
                    obj = varargin{1};
                    in_images = varargin{2};
                    method = varargin{3};
                    option = 'single';
                    relaxFactor = 0; 
                    dis_opt = [];
                    lp_solver = 'glpk';
                otherwise
                    error('Invalid number of input arguments (should be 1, 2, 3, 4, 5, or 6)');
            end
            
            if strcmp(method, 'approx-star') || strcmp(method, 'abs-dom') || contains(method, 'relax-star')
                images = obj.reach_star_multipleInputs(in_images, method, option, relaxFactor, dis_opt, lp_solver);
            elseif strcmp(method, 'approx-zono')
                images = obj.reach_zono_multipleInputs(in_images, option);
            else
                error('Unsupported reachability method');
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
            fprintf('\nParsing a tanh layer is done successfully');        
        end
        
    end
    
    
    
end

