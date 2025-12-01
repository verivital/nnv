function test_cnn()
    % TEST_CNN - Test CNN construction and MATLAB to NNV conversion
    %
    % Tests that:
    %   1. Empty NN constructor works
    %   2. CNN can be trained using MATLAB's Deep Learning Toolbox
    %   3. matlab2nnv successfully converts MATLAB network to NNV format

    %% Test 1: NN constructor
    net = NN;

    % ASSERTION 1: Empty NN constructor works
    assert(~isempty(net), 'Empty NN constructor should work');

    %% Test 2: CNN parse
    % Get digit dataset path
    digitDatasetPath = fullfile(matlabroot, 'toolbox', 'nnet', 'nndemos', ...
        'nndatasets', 'DigitDataset');

    % Check if dataset exists
    if ~isfolder(digitDatasetPath)
        warning('Digit dataset not found at %s. Skipping CNN training test.', digitDatasetPath);
        data = struct();
        data.test_skipped = true;
        data.skip_reason = 'Digit dataset not found';
        save_test_data(data, 'test_cnn', 'results', 'subdir', 'nn');
        return;
    end

    imds = imageDatastore(digitDatasetPath, ...
        'IncludeSubfolders', true, 'LabelSource', 'foldernames');

    numTrainFiles = 750;
    [imdsTrain, imdsValidation] = splitEachLabel(imds, numTrainFiles, 'randomize');

    % ASSERTION 2: Dataset loaded successfully
    assert(~isempty(imdsTrain), 'Training dataset should be loaded');
    assert(~isempty(imdsValidation), 'Validation dataset should be loaded');

    % Define CNN architecture
    layers = [
        imageInputLayer([28 28 1])

        convolution2dLayer(3, 8, 'Padding', 'same')
        batchNormalizationLayer
        reluLayer

        maxPooling2dLayer(2, 'Stride', 2)

        convolution2dLayer(3, 16, 'Padding', 'same')
        batchNormalizationLayer
        reluLayer

        maxPooling2dLayer(2, 'Stride', 2)

        convolution2dLayer(3, 32, 'Padding', 'same')
        batchNormalizationLayer
        reluLayer

        fullyConnectedLayer(10)
        softmaxLayer
        classificationLayer];

    % Training options (1 epoch for speed)
    options = trainingOptions('sgdm', ...
        'InitialLearnRate', 0.01, ...
        'MaxEpochs', 1, ...
        'Shuffle', 'every-epoch', ...
        'ValidationData', imdsValidation, ...
        'ValidationFrequency', 30, ...
        'Verbose', false);

    % Train the network
    fprintf('Training CNN (1 epoch)...\n');
    tic;
    MatlabNet = trainNetwork(imdsTrain, layers, options);
    t_train = toc;
    fprintf('Training completed in %.1f seconds.\n', t_train);

    % ASSERTION 3: Network trained successfully
    assert(~isempty(MatlabNet), 'Network should be trained successfully');

    %% Test 3: Convert to NNV format
    tic;
    nnvNet = matlab2nnv(MatlabNet);
    t_convert = toc;

    % ASSERTION 4: Conversion to NNV works
    assert(~isempty(nnvNet), 'matlab2nnv should produce non-empty result');
    assert(isa(nnvNet, 'NN'), 'Converted network should be of type NN');

    fprintf('Conversion to NNV completed in %.3f seconds.\n', t_convert);
    fprintf('NNV network has %d layers.\n', length(nnvNet.Layers));

    % Save regression data
    data = struct();
    data.t_train = t_train;
    data.t_convert = t_convert;
    data.num_layers_matlab = numel(MatlabNet.Layers);
    data.num_layers_nnv = length(nnvNet.Layers);
    data.input_size = [28 28 1];
    data.num_classes = 10;
    save_test_data(data, 'test_cnn', 'results', 'subdir', 'nn');
end
