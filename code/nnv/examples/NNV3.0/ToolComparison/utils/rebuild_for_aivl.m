function net = rebuild_for_aivl(raw)
%REBUILD_FOR_AIVL  Produce an AIVL-acceptable dlnetwork from an imported one.
%
%   net = rebuild_for_aivl(raw)
%
%   Two modes:
%
%   (1) FC-only sequential networks (ACAS Xu, RL): walks the layer list,
%       folds each ScalingLayer with unit Scale into its preceding
%       FullyConnectedLayer's Bias, strips adapter layers, and returns a
%       clean featureInputLayer -> FC -> activation -> ... dlnetwork.
%
%   (2) CNN networks (oval21, collins_rul, MNIST-ResNet): a lighter touch
%       -- AIVL accepts conv/relu/maxpool/batchnorm/addition/dropout etc.
%       directly, but rejects VerifyBatchSizeLayer and ScalingLayer-with-
%       non-trivial-affine. So for CNN inputs we just remove those two
%       layer types and keep everything else, then rebuild via dlnetwork.
%
%   Auto-selects mode by checking whether a FullyConnectedLayer exists.

    L = raw.Layers;
    fcIdx = find(arrayfun(@(l) isa(l,'nnet.cnn.layer.FullyConnectedLayer'), L), 1);
    convIdx = find(arrayfun(@(l) isa(l,'nnet.cnn.layer.Convolution2DLayer'), L), 1);

    if isempty(convIdx)
        % Pure FC: full reconstruction (ACAS Xu, RL).
        if isempty(fcIdx)
            error('rebuild_for_aivl: no FullyConnectedLayer or Convolution2DLayer in imported network');
        end
        net = rebuild_fc_only(L, L(fcIdx).InputSize);
    else
        % CNN path: surgical strip of disallowed adapters.
        net = strip_aivl_disallowed(raw);
    end
end

function net = rebuild_fc_only(L, inSize)
    newLayers = featureInputLayer(inSize, 'Name', 'in_rebuilt');
    i = 1;
    while i <= numel(L)
        cur = L(i);
        if isa(cur, 'nnet.cnn.layer.FullyConnectedLayer')
            W = cur.Weights; b = cur.Bias;
            if i+1 <= numel(L) && isa(L(i+1), 'nnet.cnn.layer.ScalingLayer') ...
                    && all(L(i+1).Scale(:) == 1)
                b = reshape(b(:) + L(i+1).Offset(:), size(b));
                i = i + 1;
            end
            newLayers(end+1, 1) = fullyConnectedLayer(cur.OutputSize, ...
                'Name', sprintf('fc_%d', numel(newLayers)), ...
                'Weights', W, 'Bias', b); %#ok<AGROW>
        elseif isa(cur, 'nnet.cnn.layer.ReLULayer')
            newLayers(end+1, 1) = reluLayer('Name', sprintf('relu_%d', numel(newLayers))); %#ok<AGROW>
        elseif isa(cur, 'nnet.cnn.layer.TanhLayer')
            newLayers(end+1, 1) = tanhLayer('Name', sprintf('tanh_%d', numel(newLayers))); %#ok<AGROW>
        elseif isa(cur, 'nnet.cnn.layer.SigmoidLayer')
            newLayers(end+1, 1) = sigmoidLayer('Name', sprintf('sig_%d', numel(newLayers))); %#ok<AGROW>
        elseif isa(cur, 'nnet.cnn.layer.ScalingLayer')
            if ~all(L(i).Offset(:) == 0) || ~all(L(i).Scale(:) == 1)
                warning('rebuild_for_aivl: skipping non-trivial standalone ScalingLayer %s', cur.Name);
            end
        end
        i = i + 1;
    end
    net = dlnetwork(newLayers);
end

function net = strip_aivl_disallowed(raw)
% For CNN inputs, AIVL accepts most layers in R2024b/R2025b (conv, relu,
% maxpool, batchnorm, addition, dropout, fully-connected, custom-output,
% flatten-into-2d). It rejects VerifyBatchSizeLayer and ScalingLayer with
% non-unit Scale or non-zero Offset. Strip those by removing the offending
% layer node and reconnecting source -> destination across the gap.
    lg = layerGraph(raw);
    layers = lg.Layers;
    for k = 1:numel(layers)
        L = layers(k);
        drop = false;
        if isa(L, 'nnet.onnx.layer.VerifyBatchSizeLayer')
            drop = true;
        elseif isa(L, 'nnet.cnn.layer.ScalingLayer')
            if ~all(L.Scale(:) == 1) || ~all(L.Offset(:) == 0)
                warning('rebuild_for_aivl: dropping non-trivial ScalingLayer %s', L.Name);
            end
            drop = true;
        end
        if ~drop, continue; end
        srcs = lg.Connections.Source(strcmp(lg.Connections.Destination, L.Name));
        dsts = lg.Connections.Destination(strcmp(lg.Connections.Source,      L.Name));
        lg = removeLayers(lg, L.Name);
        for s = 1:numel(srcs)
            for d = 1:numel(dsts)
                lg = connectLayers(lg, srcs{s}, dsts{d});
            end
        end
    end
    net = dlnetwork(lg);
end
