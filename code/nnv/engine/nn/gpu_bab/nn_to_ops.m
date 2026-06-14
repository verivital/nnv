function ops = nn_to_ops(nnvnet)
% NN_TO_OPS  Flatten an NNV NN object into an affine/relu op list for the GPU-BaB
%   engine. Phase 1 supports feedforward FullyConnected + ReLU nets (acasxu,
%   MNIST-FC). Input/Flatten layers are skipped (no effect on a flat input box);
%   any other layer ERRORS, so coverage is explicit rather than silently wrong.
%
%   ops: cell array of structs consumed by gpu_bab_ibp / the CROWN pass:
%     struct('type','affine','W',W,'b',b)   % y = W*x + b
%     struct('type','relu')                 % y = max(0,x)
    ops = {};
    for i = 1:numel(nnvnet.Layers)
        L = nnvnet.Layers{i};
        cls = class(L);
        if contains(cls, 'FullyConnected')
            ops{end+1} = struct('type', 'affine', 'W', L.Weights, 'b', L.Bias); %#ok<AGROW>
        elseif contains(cls, 'Relu') || contains(cls, 'ReLU')
            ops{end+1} = struct('type', 'relu'); %#ok<AGROW>
        elseif contains(cls, 'Input') || contains(cls, 'Flatten')
            % skip -- identity on a flat input box. NOTE: a normalizing ImageInputLayer
            % (e.g. zerocenter) is an AFFINE shift the caller must apply to the input
            % (pre-normalize the image / box) since it is dropped here.
        elseif contains(cls, 'Softmax') || contains(cls, 'Classification') ...
                || contains(cls, 'Output') || contains(cls, 'Regression')
            % skip -- softmax is monotonic and the output layer is identity on logits,
            % so neither affects argmax-robustness verification on the pre-softmax scores.
        else
            error('nn_to_ops:unsupported', ...
                ['GPU-BaB Phase 1 supports FullyConnected + ReLU only; got "%s" ' ...
                 '(layer %d). Add a backward rule for this layer type to extend coverage.'], ...
                cls, i);
        end
    end
end
