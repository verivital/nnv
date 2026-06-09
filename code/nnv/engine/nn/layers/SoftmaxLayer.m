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
            % Softmax reachability.
            %
            % Two regimes, selected by obj.IsFinalLayer:
            %  - FINAL-layer softmax (default): identity passthrough. SOUND only
            %    for argmax-style robustness specs stated on the pre-softmax
            %    logits, because softmax is per-coordinate monotone so
            %    argmax(softmax(x)) == argmax(x). NOT valid for specs on the
            %    softmax *probabilities* (e.g. P(class) > 0.9) -- use a
            %    non-final softmax (IsFinalLayer=false) for those.
            %  - INTERMEDIATE softmax (IsFinalLayer=false, e.g. inside
            %    attention): identity is UNSOUND, so compute a sound interval
            %    over-approximation (nn/funcs/Softmax.m). We deliberately DO NOT
            %    silently fall back to identity on failure -- that would hide
            %    unsoundness; an unsupported set type raises an error instead.
            %
            % NOTE: the previous implementation routed through
            % Softmax.reach_star_approx(IS, method) with method='' (no method
            % passed by most callers); that dispatch errored on the empty
            % method and a silent catch returned the UNSOUND identity. We now
            % call the bounds directly so intermediate softmax is actually sound.
            obj = varargin{1};
            IS  = varargin{2};

            if isempty(IS)
                OS = IS; return;
            end

            is_final = true;
            if isprop(obj, 'IsFinalLayer') && ~obj.IsFinalLayer
                is_final = false;
            end

            if is_final
                OS = IS;   % sound for argmax-on-logits specs (see note above)
                return;
            end

            % Intermediate softmax: sound bounds, NO silent identity fallback.
            if isa(IS, 'Star')
                OS = Softmax.reach_star_approx_bounds(IS);
            elseif isa(IS, 'ImageStar')
                S = IS.toStar();
                Sout = Softmax.reach_star_approx_bounds(S);
                OS = Sout.toImageStar(IS.height, IS.width, IS.numChannel);
            elseif isa(IS, 'Zono')
                OS = Softmax.reach_zono_approx(IS);
            else
                error('SoftmaxLayer:reach', ...
                    ['Sound intermediate-softmax reach is not implemented for set type ''%s''. ' ...
                     'Refusing to return an unsound identity passthrough.'], class(IS));
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

