classdef UpsampleLayer < handle
    % UpsampleLayer in NN
    % upsamples the input into a specific dimensions, scaled by the 'scale' dimensions
    % tipically a 2D into
    % 4D, from FullyConnected to Conv2D layers

    % author: Neelanjana Pal
    % Date: 06/20/2023
    
    properties
        Name = 'upsample_layer'; % default
        NumInputs = 1;          % default
        InputNames = {'in'};    % default
        NumOutputs = 1;         % default
        OutputNames = {'out'};  % default
        scaleDim = [1 1 1 1];   % default
    end
    
    methods
        function obj = UpsampleLayer(varargin)
            switch nargin
                
                case 6
                    name = varargin{1};
                    numInputs = varargin{2};
                    numOutputs = varargin{3};
                    inputNames = varargin{4};
                    outputNames = varargin{5};
                    scaleDim = varargin{6};
                    
                case 2
                    name = varargin{1};
                    scaleDim = varargin{2};
                    numInputs = 1;
                    numOutputs = 1;
                    inputNames = {'in'};
                    outputNames = {'out'};
                                        
                case 0
                    name = 'reshape_layer';
                    numInputs = 1;
                    numOutputs = 1;
                    inputNames = {'in'};
                    outputNames = {'out'};
                    scaleDim = [1 1 1 1];
             
                otherwise
                    
                    error('Invalid number of input arguments, should be 0, 2, or 6');        
            end
            
            
            if ~ischar(name)
                error('Invalid name, should be a charracter array');
            end
            
            if numInputs < 1
                error('Invalid number of inputs');
            end
                       
            if numOutputs < 1
                error('Invalid number of outputs');
            end
            
            if ~iscell(inputNames)
                error('Invalid input names, should be a cell');
            end
            
            if ~iscell(outputNames)
                error('Invalid output names, should be a cell');
            end
            
            obj.Name = name;
            obj.NumInputs = numInputs;
            obj.NumOutputs = numOutputs;
            obj.InputNames = inputNames;
            obj.OutputNames = outputNames;
            obj.scaleDim = flip(scaleDim);
                        
        end
    end

    methods
        
        % evaluate
        function upsample = evaluate(obj,image)
            %@image: an multi-channels image
            %@flatten_im: flatten image
            
            if length(obj.scaleDim) == 4
                scaleDim = obj.scaleDim(1:3);
            elseif length(obj.scaleDim) == 3
                scaleDim = obj.scaleDim;
            end

            upsample = dlresize(dlarray(image), 'Scale', scaleDim, 'DataFormat', "SSS");%, 'Method', "nearest", 'GeometricTransformMode',  "half_pixel", 'NearestRoundingMode', "round");
            upsample = extractdata(upsample);   
        end
    end
    
    methods % reachability method
        
        function image = reach_single_input(obj, in_image)
            % @in_image: input imagestar
            % @image: output set
            
            % TODO: implement this function, just need to modify the
            % dimensions of an ImageStar or convert a Star to an ImageStar
            % Should also support ImageZono and Zono
            
            if length(obj.scaleDim) == 4
                scaleDimC = obj.scaleDim(1:3);

            elseif length(obj.scaleDim) == 3
                scaleDimC = obj.scaleDim;
            end
            
            c = extractdata(dlresize(dlarray(in_image.V(:,:,:,1)), 'Scale', scaleDimC, 'DataFormat', "SSS"));

            if size(in_image.V,4)>2
                scaleDimV = [scaleDimC, 1];
                dataformat = "SSSS";
            else
                scaleDimV = scaleDimC;
                dataformat = "SSS";
            end

            V = extractdata(dlresize(dlarray(in_image.V(:,:,:,2:end)), 'Scale', scaleDimV, 'DataFormat', dataformat));
            V_new = cat(4,c,V);
            if isempty(in_image.im_lb) && isempty(in_image.im_ub)
                im_lb = [];
                im_ub = [];
            else
                im_lb = extractdata(dlresize(dlarray(in_image.im_lb), 'Scale', scaleDimC, 'DataFormat', "SSS"));
                im_ub = extractdata(dlresize(dlarray(in_image.im_lb), 'Scale', scaleDimC, 'DataFormat', "SSS"));
            end
            image = ImageStar(V_new,in_image.C,in_image.d,in_image.pred_lb,in_image.pred_ub,im_lb,im_ub);
        end
        
        % handle multiple inputs
        function S = reach_multipleInputs(obj, inputs, option)
            % @inputs: an array of ImageStars
            % @option: = 'parallel' or 'single'
            % @S: output ImageStar
            
            n = length(inputs);
            if isa(inputs(1), 'ImageStar')
                S(n) = ImageStar;
            elseif isa(inputs(1), 'ImageZono')
                S(n) = ImageZono;
            else
                error('Unknown input data set');
            end
          
            if strcmp(option, 'parallel')
                parfor i=1:n
                    S(i) = obj.reach_single_input(inputs(i));
                end
            elseif strcmp(option, 'single') || isempty(option)
                for i=1:n
                    S(i) = obj.reach_single_input(inputs(i));
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
            
            if strcmp(method, 'approx-star') || strcmp(method, 'exact-star') || strcmp(method, 'abs-dom') || strcmp(method, 'approx-zono') || contains(method, "relax-star")
                IS = obj.reach_multipleInputs(in_images, option);
            else
                error('Unknown reachability method');
            end
  
        end
        
    end


    methods(Static)
        % parsing method
        
        function L = parse(layer)
            % @layer: Custom Upsample Layer (created during importONNXLayers)
            % @L: constructed layer
                      
            if ~contains(class(layer), 'UpsampleLayer')
                error('Input is not a upsample layer');
            end
            
            params = layer.ONNXParams.Learnables;
            par_fields = fields(params);
            if length(par_fields) == 1
                params = struct2cell(params);
                scaleDim = extractdata(params{1});
                scaleDim = reshape(scaleDim, [1 length(scaleDim)]);

            else
                error('Parsing Reshape Layer was unsuccessful. We only support reshape layer with one Nonlearnable parameter.')
            end

            L = UpsampleLayer(layer.Name, layer.NumInputs, layer.NumOutputs, layer.InputNames, layer.OutputNames, scaleDim);
        end

    end

end

