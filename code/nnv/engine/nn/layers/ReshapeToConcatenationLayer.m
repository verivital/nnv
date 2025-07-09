classdef ReshapeToConcatenationLayer < handle
    % ReshapeToConcatenationLayer: reshapes multiple inputs and concatenates them
    %
    % Author: Navid Hashemi (updated)
    % Date: 2025-07-04
    
    properties
        Name = 'ReshapeToConcatenationLayer'
        NumInputs
        InputNames    % cell array of strings, e.g., {'in1','in2','in3'}
        NumOutputs
        OutputNames   % cell array of strings, e.g., {'output0'}
        TargetShapes  % cell array of target reshape sizes (one per input)
        ConcatDim = 3 % concatenate along channel dimension by default
    end
    
    methods
        % Constructor
        function obj = ReshapeToConcatenationLayer(name, numInputs, numOutputs, inputNames, outputNames, targetShapes, concatDim)
            obj.Name = name;
            obj.NumInputs = numInputs;
            
            % Ensure InputNames is a cell array of strings of length numInputs
            if ischar(inputNames)
                obj.InputNames = {inputNames};
            elseif iscell(inputNames)
                obj.InputNames = inputNames;
            else
                error('InputNames must be a cell array of strings or char');
            end
            if length(obj.InputNames) ~= numInputs
                error('NumInputs does not match length of InputNames');
            end
            
            obj.NumOutputs = numOutputs;
            
            % Ensure OutputNames is a cell array of strings of length numOutputs
            if nargin < 5 || isempty(outputNames)
                % Default output names if missing
                obj.OutputNames = cell(1, numOutputs);
                for i=1:numOutputs
                    obj.OutputNames{i} = sprintf('output%d', i-1); % zero-based output names like 'output0'
                end
            elseif ischar(outputNames)
                obj.OutputNames = {outputNames};
            elseif iscell(outputNames)
                obj.OutputNames = outputNames;
            else
                error('OutputNames must be a cell array of strings or char');
            end
            if length(obj.OutputNames) ~= numOutputs
                error('NumOutputs does not match length of OutputNames');
            end
            
            obj.TargetShapes = targetShapes;
            if nargin >= 7 && ~isempty(concatDim)
                obj.ConcatDim = concatDim;
            end
        end
        
        % Single input reachability (inputStars is a cell array)
        function outputStar = reach_single_input(obj, inputStars)
            if ~iscell(inputStars)
                error('reach_single_input expects a cell array of ImageStars.');
            end
            numInputs = length(inputStars);
            reshapedV = cell(1, numInputs);
            C = inputStars{1}.C;
            d = inputStars{1}.d;
            pred_lb = inputStars{1}.pred_lb;
            pred_ub = inputStars{1}.pred_ub;
            
            for k = 1:numInputs
                V = inputStars{k}.V;
                targetShape = obj.TargetShapes{k};
                numPred = size(V,4);
                reshapedSlices = cell(1, numPred);
                for i = 1:numPred
                    reshapedSlices{i} = reshape(V(:,:,:,i), targetShape);
                end
                reshapedV{k} = cat(4, reshapedSlices{:});
            end
            
            Vcat = cat(obj.ConcatDim, reshapedV{:});
            
            outputStar = ImageStar(Vcat, C, d, pred_lb, pred_ub);
        end
        
        % Reachability for multiple batches of inputs
        function S = reach_multipleInputs(obj, inputs, option)
            n = length(inputs);
            S(n) = ImageStar;
            
            if strcmp(option, 'parallel')
                parfor i = 1:n
                    S(i) = obj.reach_single_input(inputs{i});
                end
            elseif strcmp(option, 'single') || isempty(option)
                for i = 1:n
                    S(i) = obj.reach_single_input(inputs{i});
                end
            else
                error('Unknown computation option.');
            end
        end
        
        % Main reach dispatcher
        function output = reach(obj, varargin)
            narginchk(2,4);
            inputs = varargin{1};
            if length(varargin) >= 3
                option = varargin{3};
            else
                option = [];
            end
            
            if iscell(inputs) && iscell(inputs{1})
                output = obj.reach_multipleInputs(inputs, option);
            else
                output = obj.reach_single_input(inputs);
            end
        end
    end
    
    methods (Static)
        function L = parse(layer)
            if ~contains(class(layer), "Reshape_To_ConcatLayer")
                error('Layer is not a Reshape_To_ConcatLayer.');
            end
            
            % Extract target shapes from ONNXParams or fallback
            if isprop(layer, 'ONNXParams') && ~isempty(layer.ONNXParams)
                names = fieldnames(layer.ONNXParams.Nonlearnables);
                targetShapes = cell(1, length(names));
                for i = 1:length(names)
                    p = layer.ONNXParams.Nonlearnables.(names{i});
                    if isprop(p, 'Value')
                        targetShapes{i} = double(p.Value(:))';
                    else
                        targetShapes{i} = double(p(:))';
                    end
                end
            else
                % fallback to empty shapes (1 per input)
                targetShapes = repmat({[]}, 1, layer.NumInputs);
            end
            
            % Use exact InputNames and OutputNames from the original layer
            inputNames = layer.InputNames;
            outputNames = layer.OutputNames;
            if isempty(outputNames)
                % Fallback default output names
                outputNames = cell(1, layer.NumOutputs);
                for i=1:layer.NumOutputs
                    outputNames{i} = sprintf('output%d', i-1);
                end
            end
            
            L = ReshapeToConcatenationLayer(...
                layer.Name, ...
                layer.NumInputs, ...
                layer.NumOutputs, ...
                inputNames, ...
                outputNames, ...
                targetShapes, ...
                3); % default concat dimension
        end
    end
end



% classdef ReshapeToConcatLayer < handle
%     % ReshapeToConcatLayer: reshapes multiple inputs and concatenates them
%     %
%     % Author: ChatGPT
%     % Date: 2025-07-04
% 
%     properties
%         Name = 'ReshapeToConcatLayer'
%         NumInputs
%         InputNames
%         NumOutputs
%         OutputNames
%         TargetShapes  % cell array of reshape target sizes (one per input)
%         ConcatDim = 3 % concatenate along channel dimension by default
%     end
% 
%     methods
%         % Constructor
%         function obj = ReshapeToConcatLayer(name, numInputs, numOutputs, inputNames, outputNames, targetShapes, concatDim)
%             obj.Name = name;
%             obj.NumInputs = numInputs;
%             obj.NumOutputs = numOutputs;
%             obj.InputNames = inputNames;
%             obj.OutputNames = outputNames;
%             obj.TargetShapes = targetShapes;
%             if nargin >= 7 && ~isempty(concatDim)
%                 obj.ConcatDim = concatDim;
%             end
%         end
% 
%         % Single input reachability (for a set of ImageStars to be combined)
%         function outputStar = reach_single_input(obj, inputStars)
%             % inputStars: cell array of ImageStars, one per input
%             if ~iscell(inputStars)
%                 error('reach_single_input expects a cell array of ImageStars.');
%             end
%             numInputs = length(inputStars);
%             reshapedV = cell(1, numInputs);
%             C = inputStars{1}.C;
%             d = inputStars{1}.d;
%             pred_lb = inputStars{1}.pred_lb;
%             pred_ub = inputStars{1}.pred_ub;
% 
%             for k = 1:numInputs
%                 V = inputStars{k}.V;
%                 targetShape = obj.TargetShapes{k};
%                 numPred = size(V,4);
%                 reshapedSlices = cell(1, numPred);
%                 for i = 1:numPred
%                     reshapedSlices{i} = reshape(V(:,:,:,i), targetShape);
%                 end
%                 reshapedV{k} = cat(4, reshapedSlices{:});
%             end
% 
%             Vcat = cat(obj.ConcatDim, reshapedV{:});
% 
%             outputStar = ImageStar(Vcat, C, d, pred_lb, pred_ub);
%         end
% 
%         % Reachability for multiple batches of inputs
%         function S = reach_multipleInputs(obj, inputs, option)
%             % inputs: cell array of cell arrays (batch), e.g., { {in1_1, in2_1, in3_1}, {in1_2,...} }
%             n = length(inputs);
%             S(n) = ImageStar;
% 
%             if strcmp(option, 'parallel')
%                 parfor i = 1:n
%                     S(i) = obj.reach_single_input(inputs{i});
%                 end
%             elseif strcmp(option, 'single') || isempty(option)
%                 for i = 1:n
%                     S(i) = obj.reach_single_input(inputs{i});
%                 end
%             else
%                 error('Unknown computation option.');
%             end
%         end
% 
%         % Main reach dispatcher
%         function output = reach(obj, varargin)
%             % usage: reach(obj, inputStars, method, option)
%             narginchk(2,4);
%             inputs = varargin{1};
%             if length(varargin) >= 3
%                 option = varargin{3};
%             else
%                 option = [];
%             end
% 
%             if iscell(inputs) && iscell(inputs{1})
%                 % multiple batches
%                 output = obj.reach_multipleInputs(inputs, option);
%             else
%                 % single batch
%                 output = obj.reach_single_input(inputs);
%             end
%         end
%     end
% 
%     methods (Static)
%         function L = parse(layer)
%             % parse from MATLAB layer
%             if ~contains(class(layer), "Reshape_To_ConcatLayer")
%                 error('Layer is not a Reshape_To_ConcatLayer.');
%             end
% 
%             % Determine target shapes from ONNXParams or learnable parameters
%             if isprop(layer, 'ONNXParams') && ~isempty(layer.ONNXParams)
%                 names = fieldnames(layer.ONNXParams.Nonlearnables);
%                 targetShapes = cell(1, length(names));
%                 for i = 1:length(names)
%                     p = layer.ONNXParams.Nonlearnables.(names{i});
%                     if isprop(p,'Value')
%                         targetShapes{i} = double(p.Value(:))';
%                     else
%                         targetShapes{i} = double(p(:))';
%                     end
%                 end
%             else
%                 % fallback to empty shapes
%                 targetShapes = repmat({[]}, 1, layer.NumInputs);
%             end
% 
%             L = ReshapeToConcatLayer(...
%                 layer.Name, ...
%                 layer.NumInputs, ...
%                 layer.NumOutputs, ...
%                 layer.InputNames, ...
%                 layer.OutputNames, ...
%                 targetShapes, ...
%                 3); % default concat dimension
%         end
%     end
% end
