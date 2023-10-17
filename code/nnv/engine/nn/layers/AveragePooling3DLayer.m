classdef AveragePooling3DLayer < handle
    % The average pooling 3D layer class
    %   Contain constructor and reachability analysis methods
    % Main references:
    %    https://www.mathworks.com/help/deeplearning/ref/nnet.cnn.layer.averagepooling3dlayer.html
    %
    %   Diego Manzanas: 10/11/2023
    
    properties
        Name = 'average_pooling_3d_layer';
        % Hyperparameters
        PoolSize = [2 2 2];
        Stride = [1 1 1]; % step size for traversing input
        PaddingMode = 'manual'; 
        PaddingSize = [0 0 0; 0 0 0]; % size of padding 
        
        NumInputs = 1;
        InputNames = {'in'};
        NumOutputs = 1;
        OutputNames = {'out'};
        
    end
    
    
    methods % setting hyperparameters
        
        % constructor of the class
        function obj = AveragePooling3DLayer(varargin)
            % process inputs (variable)
            switch nargin
                
                case 8 % used for parsing a matlab averagePooling3DLayer
                    
                    name = varargin{1};
                    poolSize = varargin{2};
                    stride = varargin{3};
                    paddingSize = varargin{4};
                    
                    obj.NumInputs = varargin{5};
                    obj.InputNames = varargin{6};
                    obj.NumOutputs = varargin{7};
                    obj.OutputNames = varargin{8};
                    
                    if ~ischar(name)
                        error('Name is not char');
                    else
                        obj.Name = name;
                    end                    
                    
                    if size(poolSize, 1) ~= 1 || size(poolSize, 2) ~= 3
                        error('Invalid pool size');
                    else
                        obj.PoolSize = poolSize;
                    end
                    
                    if size(stride, 1) ~= 1 || size(stride, 2) ~= 3
                        error('Invalid stride');
                    else 
                        obj.Stride = stride; 
                    end
                    
                    if size(paddingSize, 1) ~= 2 || size(paddingSize, 2) ~= 3
                        error('Invalide padding size');
                    else
                        obj.PaddingSize = paddingSize;
                    end
                
                case 4
                    
                    name = varargin{1};
                    poolSize = varargin{2};
                    stride = varargin{3};
                    paddingSize = varargin{4};
                    
                    if ~ischar(name)
                        error('Name is not char');
                    else
                        obj.Name = name;
                    end                    
                    
                    if size(poolSize, 1) ~= 1 || size(poolSize, 2) ~= 3
                        error('Invalid pool size');
                    else
                        obj.PoolSize = poolSize;
                    end
                    
                    if size(stride, 1) ~= 1 || size(stride, 2) ~= 3
                        error('Invalid stride');
                    else 
                        obj.Stride = stride; 
                    end
                    
                    if size(paddingSize, 1) ~= 2 || size(paddingSize, 2) ~= 3
                        error('Invalide padding size');
                    else
                        obj.PaddingSize = paddingSize;
                    end
                
                case 3
                    
                    poolSize = varargin{1};
                    stride = varargin{2};
                    paddingSize = varargin{3};
                    
                    if size(poolSize, 1) ~= 1 || size(poolSize, 2) ~= 3
                        error('Invalid pool size');
                    else
                        obj.PoolSize = poolSize;
                    end
                    
                    if size(stride, 1) ~= 1 || size(stride, 2) ~= 3
                        error('Invalid stride');
                    else 
                        obj.Stride = stride; 
                    end
                    
                    if size(paddingSize, 1) ~= 2 || size(paddingSize, 2) ~= 3
                        error('Invalide padding size');
                    else
                        obj.PaddingSize = paddingSize;
                    end
                    
                case 0 % do ot modify default parameters
                                    
                otherwise
                    error('Invalid number of inputs (should be 0, 3 or 4)');
                                 
            end 
             
        end
        
        % set poolsize method 
        function set_poolSize(obj, poolSize)
            
            [n, m] = size(poolSize);
            
            if n ~= 1
                error('poolSize matrix shoule have one row');
            end
            
            if m == 1
                obj.PoolSize = [poolSize poolSize poolSize];
            elseif m == 3
                obj.PoolSize = poolSize;
            elseif m ~= 1 && m ~= 3
                error('Invalid poolSize matrix');
            end
            
        end
        
        % set stride method 
        function set_stride(obj, stride)
            % @stride: stride matrix
            
            [n, m] = size(stride);
            
            if n ~= 1
                error('Stride matrix shoule have one row');
            end
            
            if m == 1
                obj.Stride = [stride stride stride];
            elseif m == 3
                obj.Stride = stride;
            elseif m ~= 1 && m ~= 3
                error('Invalid stride matrix');
            end
            
        end
        
        % set padding 
        function set_padding(obj, padding)
            % @padding: padding matrix
            
            [n, m] = size(padding);
            
            if n ~= 2
                error('Padding matrix shoule have two rows');
            end
            
            if m == 1
                obj.PaddingSize = padding*ones(2,3);
            elseif m == 3
                obj.PaddingSize = padding;
            elseif m ~= 1 && m ~= 3
                error('Invalid padding matrix');
            end

        end
        
    end
        
    
    methods % evaluation functions
        
        % evaluate single input
        function y = evaluate(obj, input) 
            % @input: high-dimensional array, for example, input(:, :, :, :), 
            % @y: output
            %
            % @y: high-dimensional array (output volume), depth of output = number of filters
            
            x = dlarray(input, "SSSCB");
            y = avgpool(x, obj.PoolSize, 'Stride', obj.Stride, 'Padding', obj.PaddingSize);
            y = extractdata(y);
                   
        end

        % parallel evaluation on an array of inputs
        function y = evaluate_parallel(obj, inputs)
            % @inputs: an array of inputs
            % @y: an array of outputs 

            y = [];
            parfor i=1:length(inputs)
                y = [y obj.evaluate(inputs(i))];               
            end
            
        end
        
    end
        
    
    methods % reachability analysis functions
        
        % reachability analysis for one VolumeStar
        function S = reach_star_single_input(obj, input)
            % @inputs: an VolumeStar input set
            % @option: = 'single' single core for computation
            %          = 'parallel' multiple cores for parallel computation
            % @S: an VolumeStar with number of channels = obj.NumFilters
            
            if ~isa(input, 'VolumeStar')
                error('The input is not an VolumeStar object');
            end

            Y = obj.evaluate(input.V);                       
            S = VolumeStar(Y, input.C, input.d, input.pred_lb, input.pred_ub);

        end
        
        % handle multiple inputs
        function S = reach_star_multipleInputs(obj, inputs, option)
            % @inputs: an array of VolumeStars
            % @option: = 'parallel' or 'single'
            % @S: output VolumeStar
            
            n = length(inputs);
            S(n) = VolumeStar;

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

        % general function for reachability analysis
        function IS = reach(varargin)
            % @in_volumes: an input VolumeStar(s)
            % @IS: output set
            % @option: = 'single' or 'parallel'         
             
            switch nargin
                
                case 7
                    obj = varargin{1};
                    in_volumes = varargin{2};
                    method = varargin{3};
                    option = varargin{4};
                    % relaxFactor = varargin{5}; do not use
                    % dis_opt = varargin{6}; do not use
                    % lp_solver = varargin{7}; do not use
                
                case 6
                    obj = varargin{1};
                    in_volumes = varargin{2};
                    method = varargin{3};
                    option = varargin{4};
                    %relaxFactor = varargin{5}; do not use
                    % dis_opt = varargin{6}; do not use
                
                case 5
                    obj = varargin{1};
                    in_volumes = varargin{2};
                    method = varargin{3};
                    option = varargin{4};
                    %relaxFactor = varargin{5}; do not use
                
                case 4
                    obj = varargin{1};
                    in_volumes = varargin{2};
                    method = varargin{3};
                    option = varargin{4};
                
                case 3
                    obj = varargin{1};
                    in_volumes = varargin{2};
                    method = varargin{3};
                    option = [];
                otherwise
                    error('Invalid number of input arguments (should be 2, 3, 4, 5, or 6)');
            end
            
            if strcmp(method, 'approx-star') || strcmp(method, 'exact-star') || strcmp(method, 'abs-dom') || contains(method, "relax-star")
                IS = obj.reach_star_multipleInputs(in_volumes, option);
            elseif strcmp(method, 'approx-zono')
                IS = obj.reach_zono_multipleInputs(in_volumes, option);    
            else
                error("Unknown reachability method");
            end
            
        end
        
    end
    

    methods(Static)
        
        % parse a trained averagePooling3DLayer from matlab
        function L = parse(layer)
            % @layer: a average pooling 3d layer from matlab deep neural network tool box
            % @L : a AveragePooling3DLayer obj for reachability analysis purpose

            if ~isa(layer, 'nnet.cnn.layer.AveragePooling3DLayer')
                error('Input is not a Matlab nnet.cnn.layer.AveragePooling2DLayer class');
            end

            L = AveragePooling3DLayer(layer.Name, layer.PoolSize, layer.Stride, layer.PaddingSize, layer.NumInputs, layer.InputNames, layer.NumOutputs, layer.OutputNames);

        end
        
    end
    
    
end

