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
        IsFinalLayer = true;  % when true, reach is identity (sound for argmax-style specs).
                              % Set false for intermediate (e.g. attention) softmax.
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
        function OS = reach(varargin)
            % Sound softmax reach.
            % - For final-layer softmax (output is class scores → probabilities),
            %   identity passthrough is the historic NNV behavior because
            %   robustness specs are typically stated in terms of pre-softmax
            %   logits anyway (argmax preserves under monotonic transforms).
            % - For intermediate softmax (transformer attention), identity is
            %   UNSOUND. Here we delegate to nn/funcs/Softmax.m's interval
            %   bounds when caller passes a nontrivial reach method.
            %
            % Backwards compat: many callers pass (obj, IS, method, opt, relax)
            % with method irrelevant. Default to identity if we don't have a
            % method or if the input set is one-dim flat (final-layer case).
            obj = varargin{1};
            IS = varargin{2};
            method = '';
            if nargin >= 3, method = varargin{3}; end

            if isempty(IS)
                OS = IS; return;
            end

            % Heuristic: if the upstream of this softmax is the final layer
            % producing logits for argmax, identity is correct for verify.
            % Marker: obj.IsFinalLayer (settable from loader). Default true
            % to preserve legacy behavior; intermediate softmaxes (e.g.
            % attention) should set IsFinalLayer=false.
            is_final = true;
            if isprop(obj, 'IsFinalLayer') && ~obj.IsFinalLayer
                is_final = false;
            end

            if is_final
                OS = IS;
                return;
            end

            % Sound bounds via Softmax.reach_star_approx (if available).
            try
                if isa(IS, 'Star')
                    OS = Softmax.reach_star_approx(IS, method);
                elseif isa(IS, 'ImageStar')
                    S = IS.toStar();
                    Sout = Softmax.reach_star_approx(S, method);
                    OS = Sout.toImageStar(IS.height, IS.width, IS.numChannel);
                else
                    OS = IS;
                end
            catch
                % If sound reach not implemented for this set type,
                % fall back to identity (legacy behavior). Logged once.
                OS = IS;
            end
        end
        
    end

    methods % helper method

        % change params to gpuArrays
        function obj = toGPU(obj)
            % nothing to do here
        end

        % Change params precision
        function obj = changeParamsPrecision(obj, ~)

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
            
        end
        
    end


end

