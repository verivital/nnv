classdef LstmLayer < handle
    % The LSTMLayer layer clas
    
    %   Neelanjana Pal: 07/12/2023
    
    properties
        % Hyperparameters
        numHiddenUnits = 1;
        outputMode = []; % sequence or last
        hasStateInputs = 0;
        hasSatteOutputs = 0;
        stateActivationFunction = []; % tanh or softsign
        gateActivationFunction = []; % sigmoid or hard-sigmoid

        InputSize = 0;  % number of input
        NumInputs = 1;
        InputNames = {'in'};
        NumOutputs = 1;
        OutputNames = {'out'};  

        Name = 'lstm_layer';
        cellState = [];
        hiddenState = [];
        inputWeights = [];
        recurrentWeights = [];
        bias = [];
    end
    
    
    methods % main methods
        
        % constructor of the class
        function obj = LstmLayer(varargin)
            for i = 1 : nargin -1 
                if mod(i, 2) ~= 0
                    if strcmp(varargin{i}, 'Name')
                        obj.Name = varargin{i+1};
                    else
                        obj.Name = 'lstm_layer';
                    end
    
                    if strcmp(varargin{i},'CellState')
                        obj.cellState = varargin{i+1};
                    end
    
                    if strcmp(varargin{i},'HiddenState')
                        obj.hiddenState = varargin{i+1};
                    end
    
                    if strcmp(varargin{i},'InputWeights')
                        obj.inputWeights = varargin{i+1};
                        if isempty(obj.inputWeights)
                            error('Invalid Input Weights');
                        end
                    end
    
                    if strcmp(varargin{i},'RecurrentWeights')
                        obj.recurrentWeights = varargin{i+1};
                        if isempty(obj.recurrentWeights)
                            error('Invalid Recurrent Weights');
                        end
                    end
    
                    if strcmp(varargin{i},'Bias')
                        obj.bias = varargin{i+1};
                    end
                    
                    if strcmp(varargin{i},'NumHiddenUnits')
                        obj.numHiddenUnits = varargin{i+1};
                    end
    
                    if strcmp(varargin{i},'OutputMode')
                        obj.outputMode = varargin{i+1};
                        if isempty(obj.outputMode)
                            error('Need to specify output mode: sequence or last');
                        end
                    end
    
                    if strcmp(varargin{i},'HasStateInputs')
                        obj.hasStateInputs = varargin{i+1};
                    end
                    if strcmp(varargin{i},'HasStateOutputs')
                        obj.hasSatteOutputs = varargin{i+1};
                    end
                    if strcmp(varargin{i},'StateActivationFunction')
                        obj.stateActivationFunction = varargin{i+1};
                    end
                    if strcmp(varargin{i},'GateActivationFunction')
                        obj.gateActivationFunction = varargin{i+1};
                    end
                    if strcmp(varargin{i},'InputSize')
                        obj.InputSize = varargin{i+1};
                    end
                    if strcmp(varargin{i},'NumInputs')
                        obj.NumInputs = varargin{i+1};
                    end
                    if strcmp(varargin{i},'NumOutputs')
                        obj.NumOutputs = varargin{i+1};
                    end
                    if obj.NumInputs == 1
                        obj.InputNames = {'in'};
                    elseif obj.NumInputs == 3
                        obj.InputNames = {'in','hidden','cell'};
                    end
                    if obj.NumOutputs == 1
                        obj.OutputNames = {'out'};
                    elseif obj.NumOutputs == 3
                        obj.OutputNames = {'out','hidden','cell'};
                    end
                end
            end
            if ~ischar(obj.Name)
                obj.Name = 'lstm_layer';
            end
            if ~(size(obj.cellState, 1) == size(obj.hiddenState, 1) && size(obj.cellState,1) == obj.numHiddenUnits)
                error('Inconsistent dimension between the cell and hidden state vector');
            end
            
            if ~(size(obj.bias,1) == 4 * obj.numHiddenUnits && size(obj.bias,2) == 1)
                error('Inconsistent bias vector, dimension should be 4*numHiddenUnits X 1');
            end
            if ~(size(obj.inputWeights,1) == 4 * obj.numHiddenUnits && size(obj.inputWeights,2) == obj.InputSize)
                error('Inconsistent bias vector, dimension should be 4*numHiddenUnits X InputSize');
            end
            if ~(size(obj.recurrentWeights,1) == 4 * obj.numHiddenUnits && size(obj.recurrentWeights,2) == obj.numHiddenUnits)
                error('Inconsistent bias vector, dimension should be 4*numHiddenUnits X numHiddenUnits');
            end
        
        end
        
        % evaluation method
        function y = evaluate(obj, x)
            % @input: input
            % @y: output
            
           [m, n] = size(x); 
            if m~= obj.InputSize
                error('Inconsistent dimension of the input vector and the network input matrix Wi');
            end
            
            if n <= 0
                error('Invalid input sequence');
            end

            % initial hidden state
            h0 = dlarray(obj.hiddenState);            
            % initial cell state
            c0 = dlarray(obj.cellState);
            
            [Y,hiddenState,~] = lstm(x,h0,c0,obj.inputWeights,obj.recurrentWeights,obj.bias,"DataFormat","ST");
            if strcmp(obj.outputMode, 'sequence')
                % full hidden-state sequence (numHiddenUnits x seqLength),
                % matching MATLAB's lstmLayer OutputMode='sequence'
                % semantics and the star reachability path
                % (reach_star_single_input, outputMode == "sequence");
                % previously the first output of lstm() was discarded and
                % only the final hidden state was returned regardless of
                % OutputMode (see GitHub issue #325)
                y = extractdata(Y);
            else
                % OutputMode='last': final hidden state (numHiddenUnits x 1)
                y = extractdata(hiddenState);
            end
        end 
        
        function y = evaluateSequence(obj, input)
            % @input: input
            % @y: output
            % author: Neelanjana Pal
            % date: 1/6/2023

            y = obj.evaluate(input);
        end

        % main reachability analysis function
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
                    error('Invalid number of input arguments (should be 2, 3, 4, 5, or 6)');
            end
            
            if strcmp(method, 'exact-star')
                % the LSTM star path composes over-approximations
                % (per-step input boxing, approx-star activations,
                % McCormick product envelopes), so an exact star set is
                % not computable; fail closed rather than return an
                % over-approximation labeled exact, which NN.verify
                % would treat as a complete result
                error("LstmLayer does not support 'exact-star' reachability: the LSTM star set is an over-approximation. Use 'approx-star'.");
            elseif strcmp(method, 'approx-star') || strcmp(method, 'abs-dom') || contains(method, "relax-star")
                IS = obj.reach_star_multipleInputs(in_images, option);
            elseif strcmp(method, 'approx-zono')
                % the zonotope path (reach_zono) was non-functional and
                % had no sound semantics (see GitHub issue #326); fail
                % closed rather than risk returning a wrong set
                error('NNV:LstmLayer:zonoUnsupported', ...
                    'LstmLayer approx-zono reachability is not supported; use approx-star');
            else
                error('Unknown reachability method');
            end

        end
        
        function IS = reachSequence(varargin)
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
                    error('Invalid number of input arguments (should be 2, 3, 4, 5, or 6)');
            end
            
            if strcmp(method, 'exact-star')
                % kept accepted for backward compatibility with existing
                % scripts, but the LSTM star set composes
                % over-approximations (per-step boxing, approx-star
                % activations, McCormick product envelopes), so results
                % must not be interpreted with exact-star (complete)
                % semantics
                warning('NNV:LstmLayer:approximateReach', ...
                    "LstmLayer sequence star reachability is an over-approximation; 'exact-star' results are computed with approx-star semantics.");
            end

            if strcmp(method, 'approx-star') || strcmp(method, 'exact-star') || strcmp(method, 'abs-dom') || contains(method, "relax-star")
                IS = obj.reach_star_multipleInputs_Sequence(in_images, option);
            elseif strcmp(method, 'approx-zono')
                % the zonotope path (reach_zono) was non-functional and
                % had no sound semantics (see GitHub issue #326); fail
                % closed rather than risk returning a wrong set
                error('NNV:LstmLayer:zonoUnsupported', ...
                    'LstmLayer approx-zono reachability is not supported; use approx-star');
            else
                error('Unknown reachability method');
            end
        end
    
    end   
     
    methods % reachability methods
        
        % reachability analysis using imagestar
        function image = reach_star_single_input(obj, in_image)
            % @in_image: input imagestar
            % @image: output set
           
            
            if ~isa(in_image, 'ImageStar') && ~isa(in_image, 'Star')
                error('Input set is not an ImageStar or Star');
            end
            
            if isa(in_image, 'ImageStar')
                % reach using ImageStar (im_lb may be empty, e.g. for
                % sets produced by other layers, so check the height)
                if in_image.height ~= obj.InputSize
                    error('Inconsistency between the size of the input image and the InputSize of the network');
                end
                if ~all(obj.hiddenState == 0)
                    h_t_1 = ImageStar(obj.hiddenState,obj.hiddenState);
                else
                    h_t_1 = ImageStar();
                end
                
                if ~all(obj.cellState == 0)
                    c_t_1 = ImageStar(obj.cellState, obj.cellState);
                else
                    c_t_1 = ImageStar();
                end
                
                S = in_image.splitSet('col',in_image.width);

                for i = 1 : length(S)
                    S1 = S(i).affineMap(obj.inputWeights, obj.bias);
                    if h_t_1.isEmptySet == 0
                        S1R = h_t_1.affineMap(obj.recurrentWeights,[]);
                        S1 = S1.MinkowskiSum(S1R);
                    end
                    gates = S1.splitSet('row',4);
                    i_t = LogSig.reach(gates(1).toStar);
                    f_t = LogSig.reach(gates(2).toStar);
                    g_t = TanSig.reach(gates(3).toStar);
                    o_t = LogSig.reach(gates(4).toStar);
                    
                    c_t = i_t.HadamardProduct(g_t);
                    if c_t_1.isEmptySet == 0
                        c_t = c_t.MinkowskiSum(f_t.HadamardProduct(c_t_1.toStar));
                    end
                    ht(i) = o_t.HadamardProduct(TanSig.reach(c_t));

                    c_t_1 = c_t.toImageStar(gates(1).height,gates(1).width,gates(1).numChannel);
                    h_t_1 = ht(i).toImageStar(gates(4).height,gates(4).width,gates(4).numChannel);
                end
                if obj.outputMode == "last"
                    image = h_t_1;
                elseif obj.outputMode == "sequence"
                    lb = [];
                    ub = [];
                    for i = 1:length(ht)
                        [l, u] = ht(i).getRanges;
                        lb = [lb,l];
                        ub = [ub,u];
                    end
                    image = ImageStar(lb,ub);
                end
            end
            
        end
        
        % handle multiple inputs
        function S = reach_star_multipleInputs(obj, inputs, option)
            % @inputs: an array of ImageStars
            % @option: = 'parallel' or 'single'
            % @S: output ImageStar

            % The star-method branch of reach() dispatches here; the
            % method was previously missing, so reach() errored for
            % every star-based method.

            n = length(inputs);
            if isa(inputs, "ImageStar")
                S(n) = ImageStar;
            else
                error("Input must be ImageStar");
            end

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

        % handle multiple inputs
        function S = reach_star_multipleInputs_Sequence(obj, inputs, option)
            % @inputs: an array of ImageStars
            % @option: = 'parallel' or 'single'
            % @S: output ImageStar
            
            n = length(inputs);
            if isa(inputs, "ImageStar")
                S(n) = ImageStar;
            elseif isa(inputs, "Star")
                S(n) = Star;
            else
                error("Input must be ImageStar or Star");
            end
            
            if strcmp(option, 'parallel')
                parfor i=1:n
                    S(i) = obj.reach_star_single_input_Sequence(inputs(i));
                end
            elseif strcmp(option, 'single') || isempty(option)
                for i=1:n
                    S(i) = obj.reach_star_single_input_Sequence(inputs(i));
                end
            else
                error('Unknown computation option, should be parallel or single');
            end
            
        end
        
        function S = reach_star_single_input_Sequence(obj, inputs)
            % @inputs: input ImageStar (features x time steps)
            % @S: output set
            %
            % Symbolic star-set propagation through every time step.
            % The previous implementation evaluated the network
            % concretely: either at the im_lb/im_ub bound images, or by
            % mapping each basis vector of V through the nonlinear
            % network. Neither encloses the output of a non-monotone
            % map, so concrete trajectories sampled from inside the
            % input set routinely left the returned set.

            if ~isa(inputs, 'ImageStar')
                error('The input is not an ImageStar object');
            end
            S = obj.reach_star_single_input(inputs);
        end

        % (reachability analysis using imagezono) - NOT SUPPORTED
        % The previous reach_zono implementation was non-functional (see
        % GitHub issue #326): it referenced properties that do not exist
        % on LstmLayer (obj.Weights, obj.Bias, obj.OutputSize - the body
        % was copied from a fully-connected layer), its input-size check
        % flattened the whole sequence and compared against InputSize
        % (features per step), and even with the plumbing fixed a single
        % affine map cannot soundly enclose an LSTM. The dead code has
        % been removed and the method fails closed; a sound zonotope LSTM
        % step would need affine gate maps on zonotopes, sigmoid/tanh
        % abstractions, and a sound bilinear product enclosure.
        function image = reach_zono(obj, in_image) %#ok<STOUT,INUSD>
            error('NNV:LstmLayer:zonoUnsupported', ...
                'LstmLayer approx-zono reachability is not supported; use approx-star');
        end

        % handle mulitples inputs - NOT SUPPORTED (see GitHub issue #326)
        function S = reach_zono_multipleInputs(obj, inputs, option) %#ok<STOUT,INUSD>
            error('NNV:LstmLayer:zonoUnsupported', ...
                'LstmLayer approx-zono reachability is not supported; use approx-star');
        end
        
    end

    methods % helper method

        % change params to gpuArrays
        function obj = toGPU(obj)
            obj.bias = gpuArray(obj.bias);
            obj.inputWeights = gpuArray(obj.inputWeights);
            obj.recurrentWeights = gpuArray(obj.recurrentWeights);
        end

        % Change params precision
        function obj = changeParamsPrecision(obj, precision)
            if strcmp(precision, "double")
                obj.bias = double(obj.bias);
                obj.inputWeights = double(obj.inputWeights);
                obj.recurrentWeights = double(obj.recurrentWeights);
            elseif strcmp(precision, "single")
                obj.bias = single(obj.bias);
                obj.inputWeights = single(obj.inputWeights);
                obj.recurrentWeights = single(obj.recurrentWeights);
            else
                error("Parameter numerical precision must be 'single' or 'double'");
            end
        end
        
    end
    
    
    methods(Static)
        
        % parse a trained FullyConnectedLayer from matlab
        function L = parse(lstm_layer)
            % @lstm_layer: a lstm layer from matlab deep
            % neural network tool box
            % @L : a LSTMLayer obj for reachability analysis purpose
            
            if ~isa(lstm_layer, 'nnet.cnn.layer.LSTMLayer')
                error('Input is not a Matlab nnet.cnn.layer.LSTMLayer class');
            end
            L = LstmLayer('Name',lstm_layer.Name,'InputWeights', lstm_layer.InputWeights, ...
                'RecurrentWeights',lstm_layer.RecurrentWeights,'Bias', lstm_layer.Bias,'HiddenState', lstm_layer.HiddenState, ...
                'CellState',lstm_layer.CellState,'InputSize', lstm_layer.InputSize,'NumHiddenUnits', lstm_layer.NumHiddenUnits ...
                ,'OutputMode',lstm_layer.OutputMode,'StateActivationFunction',lstm_layer.StateActivationFunction,'GateActivationFunction', lstm_layer.GateActivationFunction, ...
               'HasStateInputs', lstm_layer.HasStateInputs,'HasStateOutputs',lstm_layer.HasStateOutputs, 'NumInputs',lstm_layer.NumInputs,'NumOutputs',lstm_layer.NumOutputs);         
            
        end
        
    end
    
end

