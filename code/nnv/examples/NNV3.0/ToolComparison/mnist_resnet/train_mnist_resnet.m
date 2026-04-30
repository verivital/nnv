function train_mnist_resnet(varargin)
%TRAIN_MNIST_RESNET Train compact ResNet classifiers for Experiment C.
%
%   Produces two dlnetwork .mat files in models/:
%     mnist_resnet8.mat  -- 28x28x1 MNIST classifier with addition layers
%     cifar_resnet8.mat  -- 32x32x3 CIFAR-10 classifier with addition layers
%
%   Both are deliberately small: enough conv depth to demonstrate residual
%   structure (so MW's additionLayer code path is exercised), small enough
%   that NNV's relax-star reach is tractable in seconds and exact-star is
%   feasible on at least MNIST.
%
%   train_resnets() trains both with default short epochs.
%   train_resnets('skipMNIST', true) skips MNIST.

    p = inputParser;
    addParameter(p, 'modelsDir',  fullfile(fileparts(mfilename('fullpath')), 'models'));
    addParameter(p, 'mnistEpochs', 3);
    addParameter(p, 'cifarEpochs', 5);
    addParameter(p, 'miniBatch',   128);
    addParameter(p, 'skipMNIST',   false);
    addParameter(p, 'skipCIFAR',   false);
    parse(p, varargin{:});
    opts = p.Results;
    if ~isfolder(opts.modelsDir), mkdir(opts.modelsDir); end

    rng(0);

    if ~opts.skipMNIST
        fprintf("\n=== Training MNIST-ResNet-8 ===\n");
        train_mnist_resnet8(opts);
    end

    if ~opts.skipCIFAR
        fprintf("\n=== Training CIFAR-ResNet-8 ===\n");
        train_cifar_resnet8(opts);
    end
end

% =========================================================================

function train_mnist_resnet8(opts)
    % Load MATLAB built-in digit dataset (28x28x1 handwritten 0-9, 10000 imgs).
    digitDatasetPath = fullfile(matlabroot,'toolbox','nnet','nndemos','nndatasets','DigitDataset');
    imds = imageDatastore(digitDatasetPath, ...
        'IncludeSubfolders', true, 'LabelSource', 'foldernames');
    [trainImds, testImds] = splitEachLabel(imds, 0.8, 'randomized');
    fprintf("  train=%d test=%d\n", numel(trainImds.Files), numel(testImds.Files));

    lgraph = mnist_resnet8_layers();
    net0   = dlnetwork(lgraph);   % R2025b's trainnet wants dlnetwork, not LayerGraph
    options = trainingOptions('adam', ...
        'MaxEpochs',         opts.mnistEpochs, ...
        'MiniBatchSize',     opts.miniBatch, ...
        'InitialLearnRate',  1e-3, ...
        'Shuffle',           'every-epoch', ...
        'Verbose',           true, ...
        'VerboseFrequency',  50, ...
        'Plots',             'none', ...
        'ExecutionEnvironment', 'cpu');

    net = trainnet(trainImds, net0, 'crossentropy', options);

    % Evaluate
    % Save dlnetwork before evaluating / collecting testset so the trained
    % weights aren't lost if downstream steps fail.
    save(fullfile(opts.modelsDir, 'mnist_resnet8.mat'), 'net', '-v7.3');
    fprintf("  wrote mnist_resnet8.mat\n");

    YPred = minibatchpredict(net, testImds);
    [~, I] = max(YPred, [], 2);
    classes = categories(trainImds.Labels);
    YPredLab = categorical(classes(I));
    acc = mean(YPredLab == testImds.Labels);
    fprintf("  MNIST-ResNet-8 test accuracy: %.3f\n", acc);

    [Xtest, Ytest] = collect_testset(testImds, 100); %#ok<ASGLU>
    save(fullfile(opts.modelsDir, 'mnist_resnet8_testset.mat'), 'Xtest', 'Ytest', '-v7.3');
    fprintf("  wrote testset (%d images)\n", numel(Ytest));
end

function train_cifar_resnet8(opts)
    % MATLAB ships CIFAR-10 since R2024a as cifar10Dataset.
    try
        cifar = cifar10Dataset;
    catch
        error("CIFAR-10 dataset not available in this MATLAB; run R2024a or later.");
    end
    trainImds = cifar.TrainingImages;
    testImds  = cifar.TestImages;
    fprintf("  train=%d test=%d\n", numel(trainImds.Files), numel(testImds.Files));

    lgraph = cifar_resnet8_layers();
    net0   = dlnetwork(lgraph);
    options = trainingOptions('adam', ...
        'MaxEpochs',         opts.cifarEpochs, ...
        'MiniBatchSize',     opts.miniBatch, ...
        'InitialLearnRate',  1e-3, ...
        'Shuffle',           'every-epoch', ...
        'Verbose',           true, ...
        'VerboseFrequency',  100, ...
        'Plots',             'none', ...
        'ExecutionEnvironment', 'cpu');

    net = trainnet(trainImds, net0, 'crossentropy', options);

    save(fullfile(opts.modelsDir, 'cifar_resnet8.mat'), 'net', '-v7.3');
    fprintf("  wrote cifar_resnet8.mat\n");

    YPred = minibatchpredict(net, testImds);
    [~, I] = max(YPred, [], 2);
    classes = categories(trainImds.Labels);
    YPredLab = categorical(classes(I));
    acc = mean(YPredLab == testImds.Labels);
    fprintf("  CIFAR-ResNet-8 test accuracy: %.3f\n", acc);

    [Xtest, Ytest] = collect_testset(testImds, 100); %#ok<ASGLU>
    save(fullfile(opts.modelsDir, 'cifar_resnet8_testset.mat'), 'Xtest', 'Ytest', '-v7.3');
    fprintf("  wrote testset (%d images)\n", numel(Ytest));
end

% =========================================================================
% Architecture definitions: both networks have additionLayer (residual).

function lgraph = mnist_resnet8_layers()
% ~8 weight layers: stem conv + 3 residual blocks + avg-pool + fc.
    in = imageInputLayer([28 28 1], 'Name','in', 'Normalization','rescale-zero-one');
    stem = [
        convolution2dLayer(3, 16, 'Padding','same', 'Name','c1')
        batchNormalizationLayer('Name','bn1')
        reluLayer('Name','r1')
    ];
    rb1 = residual_block('rb1', 16, 1);
    rb2 = residual_block('rb2', 16, 1);
    rb3 = residual_block('rb3', 16, 1);
    head = [
        globalAveragePooling2dLayer('Name','gap')
        fullyConnectedLayer(10, 'Name','fc')
        softmaxLayer('Name','sm')
    ];

    lgraph = layerGraph(in);
    lgraph = addLayers(lgraph, stem);
    lgraph = connectLayers(lgraph, 'in', 'c1');
    lgraph = chain(lgraph, 'r1', rb1);
    lgraph = chain(lgraph, [rb1{end}.Name], rb2);
    lgraph = chain(lgraph, [rb2{end}.Name], rb3);
    lgraph = addLayers(lgraph, head);
    lgraph = connectLayers(lgraph, [rb3{end}.Name], 'gap');
end

function lgraph = cifar_resnet8_layers()
    in = imageInputLayer([32 32 3], 'Name','in', 'Normalization','zerocenter');
    stem = [
        convolution2dLayer(3, 32, 'Padding','same', 'Name','c1')
        batchNormalizationLayer('Name','bn1')
        reluLayer('Name','r1')
    ];
    rb1 = residual_block('rb1', 32, 1);
    rb2 = residual_block('rb2', 32, 1);
    rb3 = residual_block('rb3', 32, 1);
    head = [
        globalAveragePooling2dLayer('Name','gap')
        fullyConnectedLayer(10, 'Name','fc')
        softmaxLayer('Name','sm')
    ];

    lgraph = layerGraph(in);
    lgraph = addLayers(lgraph, stem);
    lgraph = connectLayers(lgraph, 'in', 'c1');
    lgraph = chain(lgraph, 'r1', rb1);
    lgraph = chain(lgraph, [rb1{end}.Name], rb2);
    lgraph = chain(lgraph, [rb2{end}.Name], rb3);
    lgraph = addLayers(lgraph, head);
    lgraph = connectLayers(lgraph, [rb3{end}.Name], 'gap');
end

function blk = residual_block(name, ch, stride)
% Returns a cell array of layers forming a residual block:
%   <name>_c1 -> <name>_bn1 -> <name>_r1 -> <name>_c2 -> <name>_bn2 -> <name>_add(input)
% The skip connection is added via the calling code: connectLayers(<input>, <name>_add).
    blk = {
        convolution2dLayer(3, ch, 'Padding','same','Stride',stride, 'Name', [name '_c1'])
        batchNormalizationLayer('Name', [name '_bn1'])
        reluLayer('Name', [name '_r1'])
        convolution2dLayer(3, ch, 'Padding','same', 'Name', [name '_c2'])
        batchNormalizationLayer('Name', [name '_bn2'])
        additionLayer(2, 'Name', [name '_add'])
        reluLayer('Name', [name '_r2'])};
end

function lgraph = chain(lgraph, prevName, blk)
% Append residual block `blk` to lgraph after `prevName`. Wire the skip
% connection from prevName into the additionLayer's second input.
    blkArr = vertcat(blk{:});
    lgraph = addLayers(lgraph, blkArr);
    lgraph = connectLayers(lgraph, prevName, blk{1}.Name);
    addName = '';
    for i = 1:numel(blk)
        if isa(blk{i}, 'nnet.cnn.layer.AdditionLayer')
            addName = blk{i}.Name; break;
        end
    end
    lgraph = connectLayers(lgraph, prevName, [addName '/in2']);
end

% =========================================================================

function [Xtest, Ytest] = collect_testset(imds, n)
% Read up to `n` images into a 4-D array (HxWxCxN) + label vector.
% Promote grayscale 2-D images to 3-D HxWx1 so the (:,:,:,i) slot has rank 3.
    n = min(n, numel(imds.Files));
    files = imds.Files(1:n);
    labels = imds.Labels(1:n);
    sample = readimage(imds, 1);
    if ndims(sample) == 2, sample = reshape(sample, size(sample,1), size(sample,2), 1); end
    H = size(sample,1); W = size(sample,2); C = size(sample,3);
    Xtest = zeros(H, W, C, n, 'like', sample);
    for i = 1:n
        img = imread(files{i});
        if ndims(img) == 2, img = reshape(img, size(img,1), size(img,2), 1); end
        if C == 1 && size(img,3) == 3, img = rgb2gray(img); img = reshape(img, size(img,1), size(img,2), 1); end
        if C == 3 && size(img,3) == 1, img = repmat(img,1,1,3); end
        Xtest(:,:,:,i) = img;
    end
    Ytest = labels;
end
