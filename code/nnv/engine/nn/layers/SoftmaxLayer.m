classdef SoftmaxLayer < handle
    % Softmax Layer object
    % Author: Dung Tran
    % Date: 4/14/2020
    
    properties
        Name = 'SoftmaxLayer';
        NumInputs = 1; % default
        InputNames = {'in'}; % default
        NumOutputs = 1; % default
        OutputNames = {'out'}; % default
    end
    
    methods
        
        function obj = SoftmaxLayer(varargin)
            % @name: name of the layer
            % @NumInputs: number of inputs
            % @NumOutputs: number of outputs,
            % @InputNames: input names
            % @OutputNames: output names
            
            % author: Dung Tran
            % date:4/14/2020
            
            switch nargin
                
                case 5
                    name = varargin{1};
                    numInputs = varargin{2};
                    numOutputs = varargin{3};
                    inputNames = varargin{4};
                    outputNames = varargin{5};
                    
                case 1
                    name = varargin{1};
                    numInputs = 1;
                    numOutputs = 1;
                    inputNames = {'in'};
                    outputNames = {'out'};
                                        
                case 0
                    name = 'SoftmaxLayer';
                    numInputs = 1;
                    numOutputs = 1;
                    inputNames = {'in'};
                    outputNames = {'out'};
             
                otherwise
                    
                    error('Invalid number of input arguments, should be 0, 1, or 5');        
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
                        
        end
            
    end
        
        
    methods
        
        % evaluate
        function prob_im = evaluate(~,image)
            %@image: an multi-channels image
            % @prop: probability image
            
            % author: Dung Tran
            % date: 4/14/2020
            
            
            if isvector(image)
                prob_im = softmax(image);                
            else
                n = size(image);
                
                if n(1) == 0
                    error('Invalid input image');
                else
                    
                    if length(n) == 3
                       dlX = dlarray(image, 'SSC');
                       dlY = softmax(dlX);
                       prob_im = extractdata(dlY);
                    else
                        error('Input image is not a multi-channel image');
                    end

                end

            end
            
        end
        
        % reachability with imagestar
        function OS = reach(~, IS, ~, ~, ~)
            % @IS: imageStar input set
            % @seg_im: segmentation image
            % @method: reachability method
            % @reachOption: 'parallel' or 'single'
            % @relaxFactor: relaxation factor for reachability
            % @OS: imageStar output set = IS
            % we don't care method and reach option for softmax
            
            
            % author: Dung Tran
            % date: 4/14/2020
            
            % neglect softmax in reachability analysis 
            OS = IS; 
            
        end
        
        
        
    end
    
    
    methods(Static)
        % parsing method
        
        function L = parse(softmax_layer)
            % @softmax_layer: 
            % @L: constructed layer
                        
            % author: Dung Tran
            % date: 4/14/2020
            
            
            if ~isa(softmax_layer, 'nnet.cnn.layer.SoftmaxLayer')
                error('Input is not a Matlab nnet.cnn.layer.SoftmaxLayer');
            end
            
            L = SoftmaxLayer(softmax_layer.Name, softmax_layer.NumInputs, softmax_layer.NumOutputs, softmax_layer.InputNames, softmax_layer.OutputNames);
            fprintf('\nParsing a Matlab softmax layer is done successfully');
            
        end
        
        
    end
end

