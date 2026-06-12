function [prodNet, ok, reason] = multinet_product_net(nnvnet)
    % [prodNet, ok, reason] = multinet_product_net(nnvnet)
    %
    % Build the weight-shared PRODUCT network for the VNN-LIB 2.0 multi-network
    % `equal-to` case (both declare-network map to the SAME onnx, hence the same
    % NNV net f):
    %
    %       h([x_f; x_g]) = [f(x_f); f(x_g)]
    %
    % so a SINGLE NNV reach over the stacked input Star keeps the f/g output
    % correlation (shared predicates) that two independent reaches would lose
    % (VNNLIB2_SUPPORT_PLAN.md section 3.1).
    %
    % LAYER WHITELIST (soundness, -150 rule): only layer types whose product
    % construction is EXACTLY equivalent to running f on each half are accepted:
    %
    %   - FullyConnectedLayer (W [m x n], b [m x 1]): the stacked affine map is
    %     the block diagonal
    %         [W 0; 0 W] * [x1; x2] + [b; b] = [W*x1 + b; W*x2 + b]
    %     i.e. blkdiag(W, W) with bias [b; b] -- exact, no approximation.
    %   - ReluLayer: ReLU is ELEMENTWISE, so relu([x1; x2]) = [relu(x1); relu(x2)]
    %     and the layer passes through unchanged (a fresh ReluLayer with the same
    %     name is used to avoid sharing handle objects between two NN objects).
    %
    % Anything else (conv/pooling with spatial structure, normalization across
    % features, softmax, attention, ...) is NOT stack-equivariant in this flat
    % encoding -> ok = false, the caller must emit `unknown`.
    %
    % NOTE: the product net is assembled WITHOUT a Connections table (plain
    % sequential execution of nnvnet.Layers in order). If the source network's
    % Layers order were not its execution order, the product would diverge from
    % f-applied-twice; verify_multinet guards this with a concrete
    % evaluate-consistency cross-check before trusting any reach result.

    prodNet = []; ok = false; reason = '';

    if ~isa(nnvnet, 'NN')
        reason = 'input is not an NN object'; return;
    end
    L = nnvnet.Layers;
    if isempty(L) || ~iscell(L)
        reason = 'network has no layer cell array'; return;
    end

    pL = cell(1, numel(L));
    for i = 1:numel(L)
        Li = L{i};
        if isa(Li, 'FullyConnectedLayer')
            if ~isempty(Li.weightPerturb)
                reason = 'FullyConnectedLayer with weightPerturb is not supported in the product construction';
                return;
            end
            W = double(Li.Weights);
            b = double(Li.Bias); b = b(:);
            if isempty(W) || size(W, 1) ~= numel(b)
                reason = 'malformed FullyConnectedLayer (empty weights or bias mismatch)';
                return;
            end
            pL{i} = FullyConnectedLayer([char(Li.Name) '_x2'], blkdiag(W, W), [b; b]);
        elseif isa(Li, 'ReluLayer')
            pL{i} = ReluLayer(char(Li.Name));
        else
            reason = ['unsupported layer for the product construction: ' class(Li)];
            return;
        end
    end

    prodNet = NN(pL);
    prodNet.Name = [char(nnvnet.Name) '_product_x2'];
    ok = true;
end
