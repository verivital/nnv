classdef MaxUnpooling2DLayer < handle
    % The MaxUnPooling 2D layer class in CNN
    %   Contain constructor and reachability analysis methods
    % Main references:
    % 1) 
    %    

    
    %   Dung Tran: 4/14/2020
    
    properties
        Name = 'max_unpooling_2d_layer';
        NumInputs = 3;
        InputNames = {'in',  'indices' , 'size'};
        NumOutputs = 1;
        OutputNames = {'out'};
    end
    
    
    methods
        
        % constructor of the class
        function obj = MaxUnpooling2DLayer(varargin)           
            % author: Dung Tran
            % date: 4/14/2020    
            % update: 
            
            switch nargin
                
                case 5
                    
                    name = varargin{1};
                    numInputs = varargin{2};
                    inputNames = varargin{3};
                    numOutputs = varargin{4};
                    outputNames = varargin{5};
                    
                    if ~ischar(name)
                        error('Name is not char');
                    end                    
                    
                    if numInputs < 1 
                        error('Invalid number of inputs');
                    end
                    
                    if ~iscell(inputNames)
                        error('InputNames should be a cell');
                    end
                    
                    if numOutputs < 1
                        error('Invalid number of outputs');
                    end
                    
                    if ~iscell(outputNames)
                        error('OutputNames should be a cell');
                    end
                    
                    
                    obj.Name = name;
                    obj.NumInputs = numInputs;
                    obj.InputNames = inputNames;
                    obj.NumOutputs = numOutputs;
                    obj.OutputNames = outputNames;
                    
                case 3
                    
                    name = varargin{1};
                    numInputs = varargin{2};
                    inputNames = varargin{3};
                    
                    if ~ischar(name)
                        error('Name is not char');
                    end                    
                    
                    if numInputs < 1 
                        error('Invalid number of inputs');
                    end
                    
                    if ~iscell(inputNames)
                        error('InputNames should be a cell');
                    end
                    
                    obj.Name = name;
                    obj.NumInputs = numInputs;
                    obj.InputNames = inputNames;
                    
                case 2
                    
                    name = 'max_unpooling_2d_layer';
                    numInputs = varargin{1};
                    inputNames = varargin{2};
                    
                    if ~ischar(name)
                        error('Name is not char');
                    end                    
                    
                    if numInputs < 1 
                        error('Invalid number of inputs');
                    end
                    
                    if ~iscell(inputNames)
                        error('InputNames should be a cell');
                    end
                    
                    obj.Name = name;
                    obj.NumInputs = numInputs;
                    obj.InputNames = inputNames;
                    
                case 0
                    
                    obj.Name = 'max_unpooling_2d_layer';
                                                        
                otherwise
                    error('Invalid number of inputs (should be 0 or 3)');
                                 
            end 
             
        end
        
        
        % evaluation 
        function y = evaluate(~, input, indx, outputSize)
            % @input: input image
            % @indx: max index
            % @outputSize: inputSie
            % @y: output image
            
            % author: Dung Tran
            % date:4/19/2020
            
            dlX = dlarray(input, 'SSC');
            dlY = maxunpool(dlX, indx, outputSize); 
            y = extractdata(dlY);
            
        end
        
    end
    
    
    methods(Static)
       
        
        
        % parse a trained MaxUnPooling2dLayer from matlab
        function L = parse(max_unpooling_2d_layer)
            % @max_unpooling_2d_Layer: an MaxUnPooling2DLayer from matlab deep
            % neural network tool box
            % @L : an MaxUnPooling2DLayer obj for reachability analysis purpose
            
            % author: Dung Tran
            % date: 4/14/2020
            
            
            if ~isa(max_unpooling_2d_layer, 'nnet.cnn.layer.MaxUnpooling2DLayer')
                error('Input is not a Matlab nnet.cnn.layer.MaxUnpooling2DLayer class');
            end
            
            L = MaxUnpooling2DLayer(max_unpooling_2d_layer.Name, max_unpooling_2d_layer.NumInputs, max_unpooling_2d_layer.InputNames, max_unpooling_2d_layer.NumOutputs, max_unpooling_2d_layer.OutputNames);
            fprintf('\nParsing a Matlab max pooling 2d layer is done successfully');
            
        end
        
    end
    
    
    
end

