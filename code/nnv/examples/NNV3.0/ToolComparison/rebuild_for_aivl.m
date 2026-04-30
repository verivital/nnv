function net = rebuild_for_aivl(raw)
%REBUILD_FOR_AIVL  Produce an AIVL-acceptable dlnetwork from an imported one.
%
%   net = rebuild_for_aivl(raw)
%
%   Walks a sequentially-structured dlnetwork `raw`, folds each ScalingLayer
%   with unit Scale into its preceding FullyConnectedLayer's Bias, strips
%   adapter layers (VerifyBatchSize, ImageInput, Flatten*, RegressionOutput,
%   CustomOutputLayer), and returns a clean featureInputLayer -> FC ->
%   activation -> ... dlnetwork that AIVL's estimateNetworkOutputBounds and
%   verifyNetworkRobustness accept.
%
%   AIVL in R2025b rejects ScalingLayer and VerifyBatchSizeLayer with
%   identifier aivnv:verifyNetwork:DisallowedLayers. This rebuilder replaces
%   the CAV'23 ElementwiseAffineLayer fold (which no longer matches because
%   R2025a+ importers emit ScalingLayer instead).
%
%   Limitations: assumes pure FC + activation topology (no skip connections).
%   Correct for ACAS Xu, RL, and TLLverify benchmarks.
%   ResNet experiments use a DAG-aware variant or skip the rebuild.

    L = raw.Layers;
    % Input size: infer from the first FC layer's InputSize.
    fcIdx = find(arrayfun(@(l) isa(l,'nnet.cnn.layer.FullyConnectedLayer'), L), 1);
    if isempty(fcIdx)
        error('rebuild_for_aivl: no FullyConnectedLayer in imported network');
    end
    inSize = L(fcIdx).InputSize;

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
