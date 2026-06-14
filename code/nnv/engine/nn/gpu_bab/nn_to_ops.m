function ops = nn_to_ops(nnvnet, flattenOrder)
%   nn_to_ops(nnvnet, flattenOrder) -- flattenOrder controls how the conv output tensor
%   [H W C] maps to the first FC's input columns (the conv->FC boundary). MATLAB-native
%   dlnetworks flatten column-major ('colmajor', default); ONNX-imported nets flatten
%   row-major, which mis-orders the FC columns vs the engine's column-major view (caught by
%   the degenerate-IBP==evaluate guard). The dispatcher auto-calibrates: it tries each order
%   and keeps the one whose op-list matches net.evaluate. Sound either way (wrong -> guard
%   fails -> skip). Orders: 'colmajor' | 'chw_rowmajor' (ONNX [C H W]) | 'hwc_rowmajor'.
    if nargin < 2 || isempty(flattenOrder), flattenOrder = 'colmajor'; end
% NN_TO_OPS  Flatten an NNV NN object into an op list for the GPU-BaB engine.
%   Phase 1: FullyConnected + ReLU. Phase 2 (this version): + Conv2D and the
%   ImageInputLayer normalization (folded as a leading per-element affine op), so
%   the engine covers feedforward conv nets (plain CNNs). BatchNorm / Pool /
%   Addition(residual) are NOT yet handled and ERROR (sound-by-refusal: the
%   dispatcher falls back to the Star path), so coverage is explicit, never
%   silently wrong.
%
%   ops: cell array of structs consumed by gpu_bab_ibp / the CROWN pass:
%     struct('type','affine','W',W,'b',b)                       % y = W*x + b   (flat)
%     struct('type','relu')                                     % y = max(0,x)
%     struct('type','normaffine','scale',s,'shift',t,'shape',S) % y = s.*x + t  (per-element, [H W C])
%     struct('type','conv','W',W,'b',b,'stride',..,'pad',..,'dil',..,'inShape',..,'outShape',..)
%
%   SHAPE: bounds/coeffs are carried FLAT by the engine; conv/normaffine ops reshape
%   flat <-> [H W C (B/nSpec)] (column-major) at their boundary. nn_to_ops threads the
%   running spatial shape `curShape` = [H W C] so each conv records inShape/outShape,
%   and asserts prod(curShape)==size(firstFC.W,2) at the conv->FC (flatten) boundary.
%   The reshape order is column-major; the caller MUST verify it matches reach by
%   asserting gpu_bab_ibp(ops,center,center) ~= net.evaluate(center) (the input-order
%   soundness guard) before trusting any verdict.
    ops = {};
    curShape = [];   % [H W C] once an image/conv path is entered; [] for a flat FC net
    for i = 1:numel(nnvnet.Layers)
        L = nnvnet.Layers{i};
        cls = class(L);
        if contains(cls, 'FullyConnected')
            W = L.Weights;
            if ~isempty(curShape)
                % conv -> flatten -> FC boundary: the FC consumes the flattened conv
                % output; its input width must equal prod(curShape) or the flatten
                % order is inconsistent (caught here rather than producing wrong bounds).
                if size(W,2) ~= prod(curShape)
                    error('nn_to_ops:flattenMismatch', ...
                        ['FC input width %d ~= prod(curShape) %d ([%s]); flatten order ' ...
                         'inconsistent -- cannot map conv output to this FC soundly.'], ...
                        size(W,2), prod(curShape), num2str(curShape));
                end
                % reorder the FC weight columns from the net's flatten order to the engine's
                % column-major [H W C] view (identity for 'colmajor'); validated by the guard.
                W = i_flatten_permute(W, curShape, flattenOrder);
                curShape = [];   % flattened; back to flat-vector world
            end
            ops{end+1} = struct('type','affine','W',W,'b',L.Bias); %#ok<AGROW>
        elseif contains(cls, 'Relu') || contains(cls, 'ReLU')
            ops{end+1} = struct('type','relu'); %#ok<AGROW>
        elseif contains(cls, 'ImageInput') || contains(cls, 'Image3DInput')
            % seed the spatial shape from InputSize ([H W C]) and, if the layer
            % normalizes, fold the normalization as the FIRST op (the -150 input-mismatch
            % fix: the box must be bounded in the SAME normalized space reach feeds conv).
            curShape = double(L.InputSize(:)');
            [s, t, hasNorm] = i_input_norm(L);
            if hasNorm
                ops{end+1} = struct('type','normaffine','scale',s,'shift',t,'shape',curShape); %#ok<AGROW>
            end
        elseif contains(cls, 'Flatten')
            % identity on values; the conv->FC boundary assert (above) enforces the order.
        elseif contains(cls, 'Conv2D') || (contains(cls,'Conv') && ~contains(cls,'Transposed'))
            if isempty(curShape)
                error('nn_to_ops:noInputShape', ...
                    'Conv layer %d reached without a known input [H W C] shape (need an ImageInputLayer).', i);
            end
            op = i_conv_op(L, curShape);
            curShape = op.outShape;
            ops{end+1} = op; %#ok<AGROW>
        elseif contains(cls, 'BatchNorm')
            % BatchNorm at inference is a per-CHANNEL affine y = s_c.*x + t_c (LINEAR ->
            % exact, no relaxation): fold as a per-channel normaffine ([1 1 C] broadcast
            % over [H W C], handled by the existing normaffine IBP/CROWN). Conv-path only;
            % a post-flatten BN (curShape==[]) is refused (rare; sound-by-refusal).
            if isempty(curShape)
                error('nn_to_ops:bnFlat', ...
                    'BatchNorm on a flat (post-flatten) vector not yet supported -- refused for soundness.');
            end
            [s, t] = i_bn_fold(L, curShape);
            ops{end+1} = struct('type','normaffine','scale',s,'shift',t,'shape',curShape); %#ok<AGROW>
        elseif contains(cls, 'AvgPool') || contains(cls, 'AveragePool') || contains(cls, 'GlobalAveragePool')
            if isempty(curShape)
                error('nn_to_ops:poolFlat', ...
                    'Pooling on a flat vector has no spatial shape -- refused for soundness.');
            end
            op = i_avgpool_op(L, curShape, cls);
            curShape = op.outShape;
            ops{end+1} = op; %#ok<AGROW>
        elseif contains(cls, 'MaxPool') || contains(cls, 'MaxPooling')
            if isempty(curShape)
                error('nn_to_ops:poolFlat', ...
                    'Pooling on a flat vector has no spatial shape -- refused for soundness.');
            end
            op = i_maxpool_op(L, curShape);
            curShape = op.outShape;
            ops{end+1} = op; %#ok<AGROW>
        elseif contains(cls, 'Softmax') || contains(cls, 'Classification') ...
                || contains(cls, 'Output') || contains(cls, 'Regression') ...
                || contains(cls, 'Placeholder')
            % trailing output tail -- skippable for argmax-robustness ONLY if nothing
            % computational follows.
            for jj = i+1:numel(nnvnet.Layers)
                c2 = class(nnvnet.Layers{jj});
                if contains(c2,'FullyConnected') || contains(c2,'Relu') || contains(c2,'ReLU') ...
                        || contains(c2,'Conv') || contains(c2,'BatchNorm') || contains(c2,'Pool') ...
                        || contains(c2,'Addition')
                    error('nn_to_ops:unsupported', ...
                        'Layer "%s" (layer %d) is not a trailing output layer ("%s" follows).', cls, i, c2);
                end
            end
        else
            error('nn_to_ops:unsupported', ...
                ['GPU-BaB Phase 2 supports ImageInput(norm)+Conv2D+FullyConnected+ReLU+Flatten; ' ...
                 'got "%s" (layer %d). BatchNorm/Pool/Addition not yet handled -- refused for soundness.'], cls, i);
        end
    end
end

% ---- conv layer -> op (params copied in NNV-native [fh fw Cin Cout] form) -------------
function op = i_conv_op(L, inShape)
    W = L.Weights;
    b = L.Bias;
    if isempty(b), b = zeros(1,1,size(W,4), 'like', W); end
    stride = i_pair(L.Stride, [1 1]);
    pad    = i_quad(L.PaddingSize, [0 0 0 0]);     % [t b l r]
    dil    = i_pair(L.DilationFactor, [1 1]);
    ng = 1;
    if isprop(L,'NumGroups') && ~isempty(L.NumGroups) && isnumeric(L.NumGroups)
        ng = double(L.NumGroups);
    end
    if ng ~= 1
        error('nn_to_ops:groupedConv', 'grouped conv (NumGroups=%g) not yet supported -- refused for soundness.', ng);
    end
    Hin = inShape(1); Win = inShape(2);
    fh = size(W,1); fw = size(W,2);
    Hout = floor((Hin + pad(1) + pad(2) - ((fh-1)*dil(1) + 1))/stride(1) + 1);
    Wout = floor((Win + pad(3) + pad(4) - ((fw-1)*dil(2) + 1))/stride(2) + 1);
    Cout = size(W,4);
    op = struct('type','conv','W',W,'b',b,'stride',stride,'pad',pad,'dil',dil, ...
                'inShape',inShape, 'outShape',[Hout Wout Cout]);
end

% ---- BatchNorm (inference) -> per-channel affine fold -------------------------------
function [s, t] = i_bn_fold(L, shape)
% y = Scale.*(x - TrainedMean)./sqrt(TrainedVariance + Epsilon) + Offset = s_c.*x + t_c.
% Per channel -> [1 1 C] (broadcast over [H W C] by the normaffine op). LINEAR -> exact.
    C = shape(3);
    mu = L.TrainedMean(:); v = L.TrainedVariance(:);
    ga = L.Scale(:); be = L.Offset(:);
    if isscalar(L.Epsilon), ep = L.Epsilon * ones(C,1); else, ep = L.Epsilon(:); end
    if numel(mu)~=C || numel(v)~=C || numel(ga)~=C || numel(be)~=C
        error('nn_to_ops:bnChannels', ...
            'BatchNorm param lengths (mean %d var %d scale %d offset %d) ~= channels %d.', ...
            numel(mu), numel(v), numel(ga), numel(be), C);
    end
    den = sqrt(v + ep);
    s = reshape(ga ./ den, [1 1 C]);
    t = reshape(be - ga .* mu ./ den, [1 1 C]);
end

% ---- Average pooling -> depthwise non-overlapping LINEAR pool ------------------------
function op = i_avgpool_op(L, inShape, cls)
% Per-channel mean over each window (all +weights -> monotone IBP, exact CROWN adjoint =
% uniform distribute). NON-OVERLAPPING only (stride==pool) and UNPADDED -- padded /
% overlapping avgpool is refused (sound-by-refusal); GlobalAvgPool is the whole H x W.
    Hin = inShape(1); Win = inShape(2); C = inShape(3);
    if contains(cls, 'Global')
        pool = [Hin Win]; stride = [Hin Win]; pad = [0 0 0 0];
    else
        pool   = i_pair(L.PoolSize, [1 1]);
        stride = i_pair(L.Stride,   pool);          % default stride = pool size
        pad    = i_quad(L.PaddingSize, [0 0 0 0]);
    end
    if any(pad ~= 0)
        error('nn_to_ops:avgpoolPad', 'padded average pooling not yet supported -- refused for soundness.');
    end
    if ~isequal(stride, pool)
        error('nn_to_ops:avgpoolOverlap', 'overlapping average pooling (stride~=pool) not yet supported -- refused for soundness.');
    end
    Hout = floor((Hin - pool(1))/stride(1) + 1);
    Wout = floor((Win - pool(2))/stride(2) + 1);
    op = struct('type','avgpool','pool',pool,'stride',stride,'pad',pad, ...
                'inShape',inShape, 'outShape',[Hout Wout C]);
end

% ---- Max pooling -> depthwise per-window max (NONLINEAR; sound relax in CROWN) --------
function op = i_maxpool_op(L, inShape)
% Per-channel max over each window. Monotone (exact IBP); CROWN uses a sound relaxation
% from the window's input bounds. Overlapping (stride<pool) is allowed (the CROWN scatter
% handles it). UNPADDED only -- padded max pooling (pads with -Inf) is refused for now.
    Hin = inShape(1); Win = inShape(2); C = inShape(3);
    pool   = i_pair(L.PoolSize, [1 1]);
    stride = i_pair(L.Stride,   pool);
    pad    = i_quad(L.PaddingSize, [0 0 0 0]);
    if any(pad ~= 0)
        error('nn_to_ops:maxpoolPad', 'padded max pooling not yet supported -- refused for soundness.');
    end
    Hout = floor((Hin - pool(1))/stride(1) + 1);
    Wout = floor((Win - pool(2))/stride(2) + 1);
    op = struct('type','maxpool','pool',pool,'stride',stride,'pad',pad, ...
                'inShape',inShape, 'outShape',[Hout Wout C]);
end

% ---- conv->FC flatten: reorder FC weight columns to the engine's column-major view -----
function W2 = i_flatten_permute(W, shape, order)
% The engine flattens the conv output [H W C] COLUMN-MAJOR (H fastest). If the net flattens
% it in another order, the FC weight columns are mis-aligned. Build the permutation that maps
% each column-major position (h,w,c) to the net's flatten index, then reorder W's columns.
    H = shape(1); Wd = shape(2); C = shape(3);
    [hh, ww, cc] = ndgrid(1:H, 1:Wd, 1:C);              % linearized = column-major over [H W C]
    h = hh(:); w = ww(:); c = cc(:);
    switch order
        case 'colmajor'
            W2 = W; return;                              % already column-major -- no reorder
        case 'chw_rowmajor'                              % ONNX [C H W] flattened row-major (W fastest)
            netidx = (c-1)*H*Wd + (h-1)*Wd + w;
        case 'hwc_rowmajor'                              % [H W C] row-major (C fastest)
            netidx = (h-1)*Wd*C + (w-1)*C + c;
        otherwise
            error('nn_to_ops:flattenOrder', 'unknown flattenOrder "%s".', order);
    end
    W2 = W(:, netidx);
end

% ---- ImageInputLayer normalization -> per-element affine (scale, shift) ---------------
function [s, t, hasNorm] = i_input_norm(L)
    s = 1; t = 0; hasNorm = false;
    if ~isprop(L,'Normalization') || isempty(L.Normalization), return; end
    nrm = lower(char(string(L.Normalization)));
    switch nrm
        case {'none',''}
            return;
        case 'zerocenter'
            t = -L.Mean; s = 1; hasNorm = true;
        case 'zscore'
            sd = L.StandardDeviation; s = 1 ./ sd; t = -L.Mean ./ sd; hasNorm = true;
        case 'rescale-zero-one'
            mn = L.Min; mx = L.Max; s = 1 ./ (mx - mn); t = -mn ./ (mx - mn); hasNorm = true;
        case 'rescale-symmetric'
            mn = L.Min; mx = L.Max; s = 2 ./ (mx - mn); t = -2*mn ./ (mx - mn) - 1; hasNorm = true;
        otherwise
            error('nn_to_ops:unknownNorm', 'unhandled ImageInputLayer Normalization "%s" -- refused for soundness.', nrm);
    end
end

function v = i_pair(x, def)
    if isempty(x), v = def; elseif isscalar(x), v = [x x]; else, v = double(x(:)'); v = v(1:2); end
end
function v = i_quad(x, def)
    if isempty(x), v = def;
    elseif isscalar(x), v = [x x x x];
    elseif numel(x)==2, v = [x(1) x(1) x(2) x(2)];
    else, v = double(x(:)'); v = v(1:4); end
end
