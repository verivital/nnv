function net = load_nnv_from_mat(mat_path)
%LOAD_NNV_FROM_MAT  Build an NNV NN object from a .mat manifest produced by
% tools/onnx2nnv_python/onnx2nnv.py.  This bypasses MATLAB's
% importNetworkFromONNX entirely so we don't lose layer structure to the
% custom-layer fold-up problem.
%
% The .mat must contain:
%   layers       — cell-array (or struct array) of layer specs, each with
%                  fields: type, name, attrs, inputs, outputs, weight_keys
%   weights      — struct of named numeric arrays
%   input_name, input_shape, input_dtype, output_name, opset, src_path, checksum
%
% Returns: an NNV `NN` object with `name2indx` populated, ready for reach().

    s = load(mat_path);
    if ~isfield(s, 'layers') || ~isfield(s, 'weights')
        error('load_nnv_from_mat: %s missing layers/weights', mat_path);
    end

    layers_raw = s.layers;
    weights = s.weights;

    % Normalize layers_raw to a cell array of structs
    if iscell(layers_raw)
        layer_cells = layers_raw;
    elseif isstruct(layers_raw)
        layer_cells = num2cell(layers_raw);
    else
        error('Unexpected layers field type: %s', class(layers_raw));
    end

    n = numel(layer_cells);
    nnv_layers = cell(1, n);
    layer_names = strings(1, n);

    % First pass: instantiate each NNV layer object
    for i = 1:n
        L = layer_cells{i};
        if iscell(L), L = L{1}; end
        type_str = char(L.type);
        name_str = char(L.name);
        attrs = L.attrs;
        wkeys = normalize_keys(L.weight_keys);

        layer_names(i) = name_str;
        nnv_layers{i} = build_layer(type_str, name_str, attrs, wkeys, weights, L);
    end

    % Build connections from inputs/outputs maps
    Connections = build_connections_table(layer_cells, layer_names);

    % Construct NN object
    net = NN(nnv_layers, Connections);
    net.name2indx = containers.Map(cellstr(layer_names), 1:n);
end


function out = normalize_keys(k)
% Normalize various forms scipy.io.savemat produces for a list-of-strings:
%   - char row (one key)
%   - char matrix (one key per row, padded)
%   - cell array of char/str
%   - empty
% Returns a cellstr column.
    if isempty(k), out = {}; return; end
    if iscell(k)
        out = cellfun(@(s) strtrim(char(s)), k, 'UniformOutput', false);
        return;
    end
    if isstring(k)
        out = cellstr(strtrim(k));
        return;
    end
    if ischar(k)
        if size(k,1) > 1
            out = cell(size(k,1),1);
            for r = 1:size(k,1), out{r} = strtrim(k(r,:)); end
        else
            out = {strtrim(k)};
        end
        return;
    end
    error('Unexpected weight_keys type: %s', class(k));
end


function nnvL = build_layer(type_str, name_str, attrs, wkeys, weights, L) %#ok<INUSL>
    % Dispatch on layer type
    switch type_str
        case 'FeatureInputLayer'
            sz = attrs.InputSize;
            if isscalar(sz), in_size = double(sz); else, in_size = double(prod(sz)); end
            nnvL = FeatureInputLayer(name_str, in_size, 'none', 'auto', [], [], [], []);

        case 'ImageInputLayer'
            sz = attrs.InputSize; sz = double(sz(:)).';
            nnvL = ImageInputLayer(name_str, sz, 'none', 'auto', [], [], [], []);

        case 'FullyConnectedLayer'
            % Python side already emits W with shape [out, in] regardless
            % of original Gemm transB; bias as [out].
            W = double(weights.(wkeys{1}));
            b = double(weights.(wkeys{2}));
            b = b(:);   % column [out, 1]
            nnvL = FullyConnectedLayer(name_str, W, b);

        case 'Conv2DLayer'
            W = double(weights.(wkeys{1}));      % [F, C, kH, kW] in ONNX
            b = double(weights.(wkeys{2}));
            % NNV Conv2DLayer constructor expects MATLAB ordering [kH, kW, C, F]
            W_mat = permute(W, [3 4 2 1]);
            % NNV Conv2DLayer expects bias as [1, 1, NumFilters] — match.
            num_filters = size(W_mat, 4);
            b_mat = reshape(double(b(:)), 1, 1, num_filters);
            kshape = double(attrs.KernelSize(:)).';
            strides = double(attrs.Strides(:)).';
            pads_full = double(attrs.Pads(:)).';
            % ONNX pads: [t, l, b, r] -> MATLAB [t, l, b, r] padding rectangle
            if numel(pads_full) == 4
                pads = pads_full([1 3 2 4]);   % t b l r
            else
                pads = [0 0 0 0];
            end
            dilations = double(attrs.Dilations(:)).';
            % Constructor: Conv2DLayer(name, weights, bias, padding, stride, dilation)
            try
                nnvL = Conv2DLayer(name_str, W_mat, b_mat, pads, strides, dilations);
            catch
                nnvL = Conv2DLayer(W_mat, b_mat, pads, strides, dilations);
                nnvL.Name = name_str;
            end

        case 'BatchNormalizationLayer'
            scale = double(weights.(wkeys{1}));
            bias  = double(weights.(wkeys{2}));
            mean_ = double(weights.(wkeys{3}));
            var_  = double(weights.(wkeys{4}));
            eps_  = double(attrs.Epsilon);
            % BN.evaluate expects params shaped [1, 1, C] for HWC inputs
            nC = numel(scale);
            scale = reshape(scale(:), 1, 1, nC);
            bias  = reshape(bias(:),  1, 1, nC);
            mean_ = reshape(mean_(:), 1, 1, nC);
            var_  = reshape(var_(:),  1, 1, nC);
            try
                nnvL = BatchNormalizationLayer('Name', name_str, 'NumChannels', nC, ...
                    'TrainedMean', mean_, 'TrainedVariance', var_, ...
                    'Offset', bias, 'Scale', scale, 'Epsilon', eps_);
            catch
                nnvL = BatchNormalizationLayer(name_str, nC, mean_, var_, bias, scale, eps_);
            end

        case 'ReluLayer'
            nnvL = ReluLayer(name_str);

        case 'LeakyReluLayer'
            scale = double(attrs.Scale);
            % LeakyReluLayer constructor expects (name, NumInputs, InputNames,
            % NumOutputs, OutputNames, gamma) — 6 args.
            nnvL = LeakyReluLayer(name_str, 1, {'in'}, 1, {'out'}, scale);

        case 'SigmoidLayer'
            nnvL = SigmoidLayer(name_str);

        case 'TanhLayer'
            nnvL = TanhLayer(name_str);

        case 'SoftmaxLayer'
            nnvL = SoftmaxLayer(name_str);

        case 'ElementwiseAffineLayer'
            scale = double(weights.(wkeys{1}));
            bias_ = double(weights.(wkeys{2}));
            do_scale  = isfield(attrs,'DoScale')  && attrs.DoScale;
            do_offset = isfield(attrs,'DoOffset') && attrs.DoOffset;
            % Preserve multi-dim shape (e.g. [H,W,C] image bias for VGG-style
            % per-pixel mean subtraction). MATLAB broadcasting will handle
            % vector vs. multi-dim cases. Only flatten if it's already a
            % vector to avoid surprising shape mismatches.
            if isvector(scale) || isscalar(scale), scale = scale(:); end
            if isvector(bias_) || isscalar(bias_), bias_ = bias_(:); end
            nnvL = ElementwiseAffineLayer(name_str, scale, bias_, do_scale, do_offset);

        case 'AdditionLayer'
            % AdditionLayer(name, NumInputs, NumOutputs, InputNames, OutputNames)
            nnvL = AdditionLayer(name_str, 2, 1, {'in1','in2'}, {'out'});

        case 'FlattenLayer'
            nnvL = FlattenLayer(name_str);
            % FlattenLayer needs a Type tag for evaluate() to dispatch
            nnvL.Type = 'nnet.onnx.layer.FlattenInto2dLayer';

        case 'ReshapeLayer'
            shp = double(attrs.TargetShape(:)).';
            nnvL = ReshapeLayer(name_str, shp);
            if isfield(attrs, 'OnnxBCHW') && double(attrs.OnnxBCHW) ~= 0
                nnvL.OnnxBCHW = true;
            end

        case 'TransposeLayer'
            % NNV doesn't have a dedicated TransposeLayer; map to a
            % PlaceholderLayer that applies the permutation at evaluate
            % time. Reach() is still identity, so reach is unsound for
            % non-identity perms — but evaluate (used by xvalidate) is
            % now correct.
            perm = double(attrs.Perm(:)).';
            % ONNX perm is 0-indexed; convert to MATLAB's 1-indexed perm.
            mperm = perm + 1;
            nnvL = PlaceholderLayer(name_str, 'Transpose');
            if ~isempty(mperm) && ~isequal(mperm, 1:numel(mperm))
                nnvL.Perm = mperm;
            end

        case 'ConcatenationLayer'
            axis = double(attrs.Axis);
            in_rank = 0;
            if isfield(attrs, 'InRank'), in_rank = double(attrs.InRank); end
            % NNV (column-vector convention): the ONNX feature axis maps to
            % MATLAB's dim 1 (cat along the feature dim of an [F, 1] vector).
            % For ONNX rank-2 [B, F] axis=-1 → MATLAB Dim=1.
            % For ONNX rank-3 [B, T, F] axis=-1 → MATLAB Dim=1 (treating T
            % as a singleton in our column model — best-effort).
            if axis < 0
                if in_rank > 0
                    onnx_pos = axis + in_rank;     % 0..in_rank-1
                else
                    onnx_pos = max(0, axis + 4);
                end
            else
                onnx_pos = axis;
            end
            % Map ONNX dim index to MATLAB dim. ONNX 0 = batch (singleton in
            % NNV), ONNX 1+ = feature. NNV column vector has feature on dim 1.
            if in_rank == 2
                axis = (onnx_pos == 1) * 1 + (onnx_pos == 0) * 2;
                if axis == 0, axis = 1; end
            elseif in_rank == 4
                % BCHW: batch=0, C=1→Dim=3, H=2→Dim=1, W=3→Dim=2
                map4 = [3 3 1 2];   % index by onnx_pos+1
                axis = map4(min(max(onnx_pos+1,1), 4));
            else
                if axis < 0, axis = max(1, axis + 4); end
                if axis < 1, axis = 1; end
            end
            % NumInputs comes from the manifest's `inputs` field; ONNX
            % Concat is variadic (1..N inputs). The layer struct is `L`.
            if isstruct(L) && isfield(L, 'inputs')
                nIn = max(2, numel(L.inputs));
            else
                nIn = 2;
            end
            inNames = arrayfun(@(k) sprintf('in%d', k), 1:nIn, 'UniformOutput', false);
            nnvL = ConcatenationLayer(name_str, nIn, 1, inNames, {'out'}, axis);

        case 'MaxPooling2DLayer'
            kshape = double(attrs.KernelSize(:)).';
            strides = double(attrs.Strides(:)).';
            pads_full = double(attrs.Pads(:)).';
            if numel(pads_full) == 4, pads = pads_full([1 3 2 4]); else, pads = [0 0 0 0]; end
            try
                nnvL = MaxPooling2DLayer(name_str, kshape, strides, pads);
            catch
                nnvL = MaxPooling2DLayer(kshape, strides, pads);
                nnvL.Name = name_str;
            end

        case 'AveragePooling2DLayer'
            kshape = double(attrs.KernelSize(:)).';
            strides = double(attrs.Strides(:)).';
            pads_full = double(attrs.Pads(:)).';
            if numel(pads_full) == 4, pads = pads_full([1 3 2 4]); else, pads = [0 0 0 0]; end
            try
                nnvL = AveragePooling2DLayer(name_str, kshape, strides, pads);
            catch
                nnvL = AveragePooling2DLayer(kshape, strides, pads);
                nnvL.Name = name_str;
            end

        case 'GlobalAveragePooling2DLayer'
            try
                nnvL = GlobalAveragePooling2DLayer(name_str);
            catch
                nnvL = GlobalAveragePooling2DLayer();
                nnvL.Name = name_str;
            end

        case 'TransposedConv2DLayer'
            % ONNX ConvTranspose weights: [in_channels, out_channels, kH, kW]
            % NNV TransposedConv2DLayer expects [kH, kW, NumFilters, NumChannels]
            % where NumFilters=out, NumChannels=in.
            W = double(weights.(wkeys{1}));
            b = double(weights.(wkeys{2}));
            W_mat = permute(W, [3 4 2 1]);   % -> [kH, kW, out, in]
            num_filters = size(W_mat, 3);
            b_mat = reshape(double(b(:)), 1, 1, num_filters);
            strides = double(attrs.Strides(:)).';
            pads_full = double(attrs.Pads(:)).';
            if numel(pads_full) == 4, pads = pads_full([1 3 2 4]); else, pads = [0 0 0 0]; end
            try
                nnvL = TransposedConv2DLayer(name_str, W_mat, b_mat, pads, strides);
            catch
                nnvL = TransposedConv2DLayer(W_mat, b_mat, pads, strides);
                nnvL.Name = name_str;
            end

        case 'Resize2DLayer'
            % Best-effort: pass through as PlaceholderLayer if NNV's
            % Resize2DLayer constructor fails (signature varies)
            try
                nnvL = Resize2DLayer(name_str);
            catch
                nnvL = PlaceholderLayer(name_str, 'Resize');
            end

        case 'SignLayer'
            % NNV's SignLayer uses Sign.evaluate / Sign.reach machinery.
            % Constructor signature: SignLayer(gamma, mode).
            mode = 'polar_zero_to_pos_one';
            if isfield(attrs, 'Mode')
                mode = char(string(attrs.Mode));
            end
            nnvL = SignLayer(0, mode);
            nnvL.Name = name_str;

        case 'PlaceholderLayer'
            % If the placeholder represents an element-wise op we can
            % evaluate (Sign/Abs/etc.), tag it so PlaceholderLayer.evaluate
            % can dispatch on it. Otherwise it's a true no-op.
            tag = 'Identity';
            if isfield(attrs, 'OriginalOp')
                op = char(string(attrs.OriginalOp));
                if any(strcmp(op, {'Sign','Abs'}))
                    tag = op;
                end
            end
            nnvL = PlaceholderLayer(name_str, tag);

        otherwise
            error('Unsupported layer type in manifest: %s', type_str);
    end
end


function T = build_connections_table(layer_cells, layer_names)
% Build NNV-style connections table from layer specs. Handles the various
% shapes scipy.io.savemat produces for a list-of-[name,idx] pairs.
    Sources = strings(0,1);
    Destinations = strings(0,1);

    for i = 1:numel(layer_cells)
        L = layer_cells{i}; if iscell(L), L = L{1}; end
        ins_raw = L.inputs;
        if isempty(ins_raw), continue; end

        % Normalize inputs to a cell array of {producer_name, idx} pairs
        ins_pairs = normalize_inputs(ins_raw);

        type_str = char(L.type);
        is_multi = ismember(type_str, {'AdditionLayer','ConcatenationLayer'});

        for k = 1:numel(ins_pairs)
            src_name = ins_pairs{k}{1};
            if is_multi
                dest_name = sprintf('%s/in%d', char(L.name), k);
            else
                dest_name = char(L.name);
            end
            Sources(end+1,1) = string(src_name); %#ok<AGROW>
            Destinations(end+1,1) = string(dest_name); %#ok<AGROW>
        end
    end

    T = table(Sources, Destinations, 'VariableNames', {'Source','Destination'});
end

function out = normalize_inputs(ins)
% Returns cell array of {name_char, idx_double} pairs. The Python emitter
% writes inputs as a flat list of producer-name strings, with idx implied
% to be 0. scipy.io.savemat hands us those as a cell row of chars (or, for
% one-element lists, a single char).
    out = {};
    if iscell(ins)
        for r = 1:numel(ins)
            e = ins{r};
            if iscell(e) && numel(e) >= 1
                out{end+1} = {char(e{1}), 0}; %#ok<AGROW>
            elseif ischar(e) || isstring(e)
                out{end+1} = {char(e), 0}; %#ok<AGROW>
            else
                out{end+1} = {char(string(e)), 0}; %#ok<AGROW>
            end
        end
        return;
    end
    if ischar(ins)
        out{end+1} = {strtrim(ins), 0};
        return;
    end
    if isstring(ins)
        for r = 1:numel(ins), out{end+1} = {char(ins(r)), 0}; end %#ok<AGROW>
        return;
    end
    error('Unexpected inputs field type: %s', class(ins));
end
