classdef PlaceholderLayer < handle
    % Placeholder Layer object
    % This layer serves as a placeholder for all the layers that to not
    % have any effect on the verification algorithm such as Dropout or
    % RegressionLayer from MATLAB
    %
    % Author: Diego Manzanas Lopez
    % Date: 03/15/2023
    
    properties
        Name = 'NoOpLayer';
        Type = ''; % e.g. dropout
        Perm = [];  % optional MATLAB-style permutation order (1-indexed)
                    % set when this placeholder represents a Transpose op
        Constant = [];  % optional constant tensor; if non-empty, evaluate
                        % returns this value regardless of input. Used for
                        % ONNX initializer tensors fed as data inputs (e.g.
                        % CLS token to a Concat).
    end

    methods % constructor

        % create layer
        function obj = PlaceholderLayer(name, Ltype)
            % @name: name of the layer
            % @type: original layer type from MATLAB
            obj.Name = name;
            obj.Type = Ltype;
        end

    end

    methods % main methods

        % evaluate
        function out_im = evaluate(obj, inputs)
            % An ACTIVE-but-UNSUPPORTED op (tagged 'UnsupportedOp:<op>' by the
            % loader: Expand/Where/ScatterND/ArgMax/Gather/dynamic-Reshape/...)
            % transforms its input; returning the input unchanged silently
            % computes a WRONG network. Refuse loudly instead.
            if startsWith(obj.Type, 'UnsupportedOp:')
                error('PlaceholderLayer:unsupportedOp', ...
                    ['Layer ''%s'' wraps unsupported ONNX op ''%s'' -- identity ' ...
                     'evaluation would compute a different network. Refusing.'], ...
                    obj.Name, extractAfter(obj.Type, 'UnsupportedOp:'));
            end
            % If this is a constant-producing placeholder (e.g. a CLS-token
            % initializer fed to Concat), return the stored value.
            if ~isempty(obj.Constant)
                out_im = obj.Constant;
                return;
            end
            % return output = input, optionally permuted or with an
            % element-wise function applied (Sign/Abs/etc.).
            if ~isempty(obj.Perm)
                % MATLAB's permute requires perm to cover at least ndims(inputs).
                % If the perm is shorter (e.g. ONNX rank-2 perm but MATLAB
                % already promoted the tensor to rank-3), pad with identity.
                p = obj.Perm(:).';
                k = max(numel(p), ndims(inputs));
                if numel(p) < k
                    p = [p, (numel(p)+1):k];
                end
                out_im = permute(inputs, p);
                return;
            end
            % Element-wise op stored as the Type tag (set by the loader for
            % ops like Sign / Abs that NNV doesn't have a dedicated layer for).
            switch obj.Type
                case 'Sign'
                    out_im = sign(inputs);
                case 'Abs'
                    out_im = abs(inputs);
                case 'Floor'
                    out_im = floor(inputs);
                case 'Ceil'
                    out_im = ceil(inputs);
                case 'Round'
                    out_im = round(inputs);
                case 'Sin'
                    out_im = sin(inputs);
                case 'Cos'
                    out_im = cos(inputs);
                case 'Tan'
                    out_im = tan(inputs);
                case 'Exp'
                    out_im = exp(inputs);
                case 'Log'
                    out_im = log(inputs);
                case 'Sqrt'
                    out_im = sqrt(inputs);
                otherwise
                    out_im = inputs;
            end
        end

        function out_sq = evaluateSequence(~, inputs)
            out_sq = inputs;
        end
        
        % reachability analysis with multiple inputs
        function IS = reach(varargin)
            % For INERT placeholders (dropout, output layers, final softmax,
            % ...) identity is correct. For ACTIVE placeholders -- ones whose
            % evaluate() actually transforms the input (Constant, Perm/
            % Transpose, or an element-wise op Type) -- identity reach is
            % UNSOUND: evaluate and reach would disagree, corrupting every
            % downstream verdict. Compute a sound over-approximation where
            % possible, refuse loudly where not.
            obj = varargin{1};
            in_images = varargin{2};
            if ~obj.isActiveOp()
                IS = in_images;   % inert: identity is correct
                return;
            end
            IS = obj.reach_active(in_images);
        end

        % reachability analysis with multiple inputs
        function IS = reachSequence(varargin)
            % Same soundness rule as reach (sequence path).
            obj = varargin{1};
            in_seqs = varargin{2};
            if ~obj.isActiveOp()
                IS = in_seqs;
                return;
            end
            IS = obj.reach_active(in_seqs);
        end

        function tf = isActiveOp(obj)
            % True when evaluate() is NOT the identity: constant-producing,
            % permuting, an element-wise op handled in evaluate's switch, or
            % an active-but-unsupported op tagged by the loader.
            activeTypes = {'Sign','Abs','Floor','Ceil','Round','Sin','Cos', ...
                           'Tan','Exp','Log','Sqrt'};
            nonIdentityPerm = ~isempty(obj.Perm) && ~isequal(obj.Perm(:).', 1:numel(obj.Perm));
            tf = ~isempty(obj.Constant) || nonIdentityPerm || any(strcmp(obj.Type, activeTypes)) ...
                 || startsWith(obj.Type, 'UnsupportedOp:');
        end

        function OS = reach_active(obj, IS)
            % Sound reach for an ACTIVE placeholder.
            %  - Constant: the output is the stored tensor regardless of input
            %    -> the EXACT reachable set is the point {c}.
            %  - Monotone elementwise ops: per-coordinate box [f(lb), f(ub)]
            %    is a sound (and box-exact) over-approximation.
            %  - Abs: standard interval abs. Sin/Cos: [-1,1] (sound, loose).
            %  - Tan (unbounded across asymptotes) and Transpose (permuting a
            %    flattened set needs the tensor shape, which the layer does
            %    not have): REFUSE rather than return a wrong set.
            if startsWith(obj.Type, 'UnsupportedOp:')
                error('PlaceholderLayer:unsupportedOp', ...
                    ['Layer ''%s'' wraps unsupported ONNX op ''%s'' -- no sound ' ...
                     'reachability exists for it. Refusing to return an unsound set.'], ...
                    obj.Name, extractAfter(obj.Type, 'UnsupportedOp:'));
            end
            if ~isempty(obj.Constant)
                c = double(obj.Constant(:));
                OS = Star(c, c);   % exact point set
                return;
            end
            if ~isempty(obj.Perm) && ~isequal(obj.Perm(:).', 1:numel(obj.Perm))
                error('PlaceholderLayer:transposeReachNeedsShape', ...
                    ['Transpose placeholder ''%s'': permuting a flattened set ' ...
                     'requires the tensor shape, which this layer does not carry. ' ...
                     'Refusing to return an (unsound) identity set.'], obj.Name);
            end
            if isa(IS, 'Star') || isa(IS, 'ImageStar')
                [lb, ub] = IS.getRanges();
            else
                error('PlaceholderLayer:reachUnsupportedSet', ...
                    'Active op ''%s'' reach not implemented for set type %s.', obj.Type, class(IS));
            end
            lb = double(lb(:)); ub = double(ub(:));
            switch obj.Type
                case {'Floor','Ceil','Round','Exp','Sign'}   % monotone non-decreasing
                    f = str2func(lower(obj.Type));
                    olb = f(lb); oub = f(ub);
                case 'Sqrt'
                    if any(lb < 0)
                        error('PlaceholderLayer:domain', ...
                            'Sqrt placeholder ''%s'': input lower bound < 0.', obj.Name);
                    end
                    olb = sqrt(lb); oub = sqrt(ub);
                case 'Log'
                    if any(lb <= 0)
                        error('PlaceholderLayer:domain', ...
                            'Log placeholder ''%s'': input lower bound <= 0.', obj.Name);
                    end
                    olb = log(lb); oub = log(ub);
                case 'Abs'
                    olb = max(0, max(lb, -ub));     % 0 iff the interval straddles 0
                    oub = max(abs(lb), abs(ub));
                case {'Sin','Cos'}
                    % Always-sound (loose) bound: the codomain [-1, 1].
                    olb = -ones(size(lb)); oub = ones(size(ub));
                case 'Tan'
                    error('PlaceholderLayer:tanUnbounded', ...
                        ['Tan placeholder ''%s'': tan is unbounded across its ' ...
                         'asymptotes; no finite sound box exists in general. ' ...
                         'Refusing to return an unsound set.'], obj.Name);
                otherwise
                    error('PlaceholderLayer:activeOpNoReach', ...
                        'Active op ''%s'' has no sound reach implementation.', obj.Type);
            end
            if isa(IS, 'ImageStar')
                OS = ImageStar(reshape(olb, IS.height, IS.width, IS.numChannel), ...
                               reshape(oub, IS.height, IS.width, IS.numChannel));
            else
                OS = Star(olb, oub);
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
            % create NNV layer 
            L = PlaceholderLayer(layer.Name, class(layer));
        end

    end

end



