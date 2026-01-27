classdef ConcatenationLayer < handle
    % Concatenation Layer object
    % Concatenate arrays along a specific dimension
    % Typically use to concatenate 1D or 2D features  (vector/matrices, rather than images such as in DepthConcatenation)
    % see: https://www.mathworks.com/help/deeplearning/ref/nnet.cnn.layer.concatenationlayer.html
    % 
    % Author: Diego Manzanas Lopez
    % Date: 03/15/2023
    
    properties
        Name = 'ConcatLayer';
        NumInputs = 1; % default
        InputNames = {'in'}; % default
        NumOutputs = 1; % default
        OutputNames = {'out'}; % default
        Dim = 0 % unspecified, must be a positive integer
    end
    
    methods % constructor
        
        % create layer
        function obj = ConcatenationLayer(varargin)
            % @name: name of the layer
            % @NumInputs: number of inputs
            % @NumOutputs: number of outputs,
            % @InputNames: input names
            % @OutputNames: output names
            % @Dim: dimension for concatenating inputs
            
            switch nargin
                
                case 6
                    name = varargin{1};
                    numInputs = varargin{2};
                    numOutputs = varargin{3};
                    inputNames = varargin{4};
                    outputNames = varargin{5};
                    dim = varargin{6};
                case 2
                    name = varargin{1};
                    dim = varargin{2};
                    numInputs = 1;
                    numOutputs = 1;
                    inputNames = {'in'};
                    outputNames = {'out'};
                otherwise
                    error('Invalid number of input arguments, should be 2 or 6');        
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

            if ~isscalar(dim) || dim < 1
                error('Concat dimension must be a positive integer')
            end
            
            obj.Name = name;
            obj.NumInputs = numInputs;
            obj.NumOutputs = numOutputs;
            obj.InputNames = inputNames;
            obj.OutputNames = outputNames; 
            obj.Dim = dim;
        end
            
    end
        
        
    methods % main methods
        
        % evaluate
        function outputs = evaluate(obj, inputs)
            % concatenation layers takes usually two inputs, but allows many (N)
            %
            % first input is obj, the second is a cell array containing the inputs to concatenate
            % Initialize outputs as the first one
            outputs = inputs{1};
            % Concatenate the inputs 
            for k=2:length(inputs)
                outputs = cat(obj.Dim, outputs, inputs{k});
            end
                
        end
 
        % reach
        function outputs = reach_single_input(obj, inputs)
            % @inputs: cell array of ImageStar sets to concatenate
            % @outputs: concatenated ImageStar with properly combined constraints
            %
            % This method properly handles concatenation of ImageStar sets with
            % independent predicate variables by using block-diagonal constraint
            % matrices and preserving all constraints from all inputs.
            %
            % Fixed: 2024 - Previously only used constraints from one input,
            % causing soundness issues. Now properly combines all constraints.

            n = length(inputs);

            if n == 1
                % Single input - just return it
                outputs = inputs{1};
                return;
            end

            % Collect information about all inputs
            nVars = zeros(n, 1);  % number of predicate variables per input
            for i = 1:n
                nVars(i) = inputs{i}.numPred;
            end
            totalVars = sum(nVars);

            % Get spatial dimensions from first input
            % V has shape [H, W, C, nVar+1]
            sz1 = size(inputs{1}.V);

            % Compute output dimensions after concatenation
            totalDimSize = 0;
            for i = 1:n
                szI = size(inputs{i}.V);
                if obj.Dim <= 3
                    totalDimSize = totalDimSize + szI(obj.Dim);
                else
                    error('Concatenation dimension must be 1, 2, or 3');
                end
            end

            % Determine output spatial size
            outSize = sz1(1:3);
            outSize(obj.Dim) = totalDimSize;

            % Build the combined V matrix
            % new_V has shape [outH, outW, outC, totalVars+1]
            new_V = zeros([outSize, totalVars + 1], 'like', inputs{1}.V);

            % Track position along concatenation dimension and predicate index
            dimOffset = 0;
            varOffset = 0;

            for i = 1:n
                szI = size(inputs{i}.V);
                dimSizeI = szI(obj.Dim);
                nVarI = nVars(i);

                % Build index ranges for this input
                switch obj.Dim
                    case 1
                        idxH = dimOffset + (1:dimSizeI);
                        idxW = 1:szI(2);
                        idxC = 1:szI(3);
                    case 2
                        idxH = 1:szI(1);
                        idxW = dimOffset + (1:dimSizeI);
                        idxC = 1:szI(3);
                    case 3
                        idxH = 1:szI(1);
                        idxW = 1:szI(2);
                        idxC = dimOffset + (1:dimSizeI);
                end

                % Copy center (column 1 of V)
                new_V(idxH, idxW, idxC, 1) = inputs{i}.V(:,:,:,1);

                % Copy basis vectors to the correct predicate variable positions
                % Input i's predicates go to positions (varOffset+1) to (varOffset+nVarI)
                % In new_V, that's columns (varOffset+2) to (varOffset+nVarI+1)
                if nVarI > 0
                    new_V(idxH, idxW, idxC, varOffset + 2 : varOffset + nVarI + 1) = inputs{i}.V(:,:,:,2:end);
                end

                dimOffset = dimOffset + dimSizeI;
                varOffset = varOffset + nVarI;
            end

            % Build block-diagonal constraint matrix and concatenate d, bounds
            % Use blkdiag for constraints
            new_C = [];
            new_d = [];
            new_pred_lb = [];
            new_pred_ub = [];

            for i = 1:n
                if ~isempty(inputs{i}.C)
                    if isempty(new_C)
                        new_C = inputs{i}.C;
                    else
                        new_C = blkdiag(new_C, inputs{i}.C);
                    end
                    new_d = [new_d; inputs{i}.d];
                end
                new_pred_lb = [new_pred_lb; inputs{i}.pred_lb];
                new_pred_ub = [new_pred_ub; inputs{i}.pred_ub];
            end

            % Handle empty constraints case
            if isempty(new_C)
                new_C = zeros(0, totalVars);
                new_d = zeros(0, 1);
            end

            % Create output ImageStar with combined constraints
            outputs = ImageStar(new_V, new_C, new_d, new_pred_lb, new_pred_ub);

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
%                 IS = obj.reach_multipleInputs(in_images, option);
                IS = obj.reach_single_input(in_images);
            else
                error('Unknown reachability method');
            end
  
        end
        
    end

    
    methods % helper method

        % change params to gpuArrays
        function obj = toGPU(obj)
            % nothing to change in here (no params)
        end

        % Change params precision
        function obj = changeParamsPrecision(obj, ~)
            % nothing to change in here (no params)
        end
        
    end
    
    
    methods(Static)
        
        % parsing method
        function L = parse(layer)
            % create NNV layer from matlab
                      
            if ~isa(layer, 'nnet.cnn.layer.ConcatenationLayer')
                error('Input is not a concatenation layer');
            end
            L = ConcatenationLayer(layer.Name, layer.NumInputs, layer.NumOutputs, layer.InputNames, layer.OutputNames, layer.Dim);
        end

    end
end

