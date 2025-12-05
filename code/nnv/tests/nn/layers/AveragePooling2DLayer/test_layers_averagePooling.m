function test_layers_averagePooling()
    % TEST_LAYERS_AVERAGEPOOLING - Test AveragePooling2DLayer operations
    %
    % Tests that:
    %   1. compute_averageMap produces correct values
    %   2. Constructor works with various parameters
    %   3. evaluate produces correct output
    %   4. Single precision evaluation works
    %   5. Zero padding input is handled correctly
    %   6. reach_star_single_input works with ImageStar
    %   7. reach_star_single_input works with attacked image

    %% Test 1: AveragePooling2DLayer Compute averageMap
    I = [1 0 2 3; 4 6 6 8; 3 1 1 0; 1 2 2 4]; % input
    L = AveragePooling2DLayer([2 2], [2 2], [0 0 0 0]);
    averageMap = L.compute_averageMap(I);

    checker = [(I(1,1)+I(1,2)+I(2,1)+I(2,2))/4, (I(1,3)+I(1,4)+I(2,3)+I(2,4))/4; ...
               (I(3,1)+I(3,2)+I(4,1)+I(4,2))/4, (I(3,3)+I(3,4)+I(4,3)+I(4,4))/4];

    % ASSERTION 1: averageMap is correct
    assert(isequal(checker, averageMap), 'AverageMap should match expected values');

    %% Test 2: AveragePooling2DLayer constructor
    L1 = AveragePooling2DLayer('test_average_pooling_2d_layer', [2 2], [1 1], [0 0 0 0]);
    L2 = AveragePooling2DLayer();
    L3 = AveragePooling2DLayer([3 3], [1 1], [0 0 0 0]);

    % ASSERTION 2: Constructors work
    assert(~isempty(L1), 'Named constructor should work');
    assert(~isempty(L2), 'Default constructor should work');
    assert(~isempty(L3), 'Parameterized constructor should work');

    %% Test 3: AveragePooling2DLayer evaluation
    % original input volume: color image with 3 channels
    inputVol(:, :, 1) = [0 0 2 0 0; 1 2 0 2 0; 0 0 2 2 0; 0 2 2 2 2; 2 2 2 1 1];
    inputVol(:, :, 2) = [1 2 2 1 2; 2 1 2 0 2; 2 2 2 0 1; 1 1 1 0 0; 1 0 2 2 1];
    inputVol(:, :, 3) = [0 0 2 2 1; 0 2 1 1 2; 0 2 0 0 1; 0 2 1 0 1; 1 2 1 0 0];

    L = AveragePooling2DLayer([3 3], [2 2], [0 0 0 0]);
    y = L.evaluate(inputVol);

    % ASSERTION 3: Evaluate produces valid output
    assert(~isempty(y), 'Evaluate should produce non-empty output');

    for i = 1:3
        for j = 1:2
            for k = 1:2
                expected = sum(inputVol(2*j-1:2*j+1, 2*k-1:2*k+1, i), 'all') / 9;
                assert(abs(y(j, k, i) - expected) <= 10*eps, ...
                    sprintf('Evaluate output mismatch at (%d,%d,%d)', j, k, i));
            end
        end
    end

    %% Test 4: AveragePooling2DLayer evaluation - single precision
    inputVol_single = single(inputVol);
    L = AveragePooling2DLayer([3 3], [2 2], [0 0 0 0]);
    y_single = L.evaluate(inputVol_single);

    % ASSERTION 4: Single precision works
    assert(isa(y_single, 'single'), 'Single precision input should produce single precision output');

    %% Test 5: AveragePooling2DLayer get zero padding input
    inputVol2(:, :, 1) = [2 0 1 2 1; 1 0 2 2 2; 1 2 2 0 2; 1 2 0 0 1; 1 0 1 1 2];
    inputVol2(:, :, 2) = [0 0 1 0 1; 0 0 2 1 1; 1 1 0 1 1; 1 1 0 2 2; 2 1 2 0 0];
    inputVol2(:, :, 3) = [1 2 2 1 0; 2 0 0 2 0; 0 0 1 0 1; 1 2 0 2 0; 1 0 2 1 0];

    paddingSize = [1 1 1 1];
    L = AveragePooling2DLayer();
    L.set_padding(paddingSize);
    I_padded = L.get_zero_padding_input(inputVol2);

    % ASSERTION 5: Zero padding is correct
    assert(isequal(I_padded(2:end-1, 2:end-1, :), inputVol2), ...
        'Interior of padded input should match original');

    corners = [I_padded(1, 1, :), I_padded(1, end, :); I_padded(end, 1, :), I_padded(end, end, :)];
    top_bot = [I_padded(1, 2:end-1, :); I_padded(end, 2:end-1, :)];
    left_right = [I_padded(2:end-1, 1, :), I_padded(2:end-1, end, :)];

    assert(isempty(find(corners, 1)), 'Corners should be zero');
    assert(isempty(find(top_bot, 1)), 'Top/bottom padding should be zero');
    assert(isempty(find(left_right, 1)), 'Left/right padding should be zero');

    %% Test 6: AveragePooling2DLayer reach star digit one example
    load one_image.mat
    Center = one_image;
    Basis = rand(28, 28);

    V(:,:,1,1) = Center;
    V(:,:,1,2) = Basis;

    Constr_mat = [1; -1];
    Constr_vec = [1; 1];
    pred_lb = -1;
    pred_ub = 1;

    input_image = ImageStar(V, Constr_mat, Constr_vec, pred_lb, pred_ub);
    L = AveragePooling2DLayer([6 4], [4 4], [1 1 0 0]);
    output_image = L.reach_star_single_input(input_image);
    sampled_images = output_image.sample(2);

    % ASSERTION 6: Reach produces valid output
    assert(~isempty(output_image), 'reach_star_single_input should produce non-empty result');
    assert(~isempty(sampled_images), 'Sampling should produce non-empty result');

    fig1 = figure;
    subplot(1, 3, 1);
    imshow(input_image.V(:,:,1));
    title('28x28 input image');
    subplot(1, 3, 2);
    imshow(sampled_images{1, 1});
    title('7x7 1st output image');
    subplot(1, 3, 3);
    imshow(sampled_images{1, 2});
    title('7x7 2nd output image');

    save_test_figure(fig1, 'test_layers_averagePooling', 'digit_reach', 1, 'subdir', 'nn/avgpool');

    %% Test 7: AveragePooling2DLayer reach star exact single input
    IM(:, :, 1) = [0 0 2 0 0; 1 2 0 2 0; 0 0 2 2 0; 0 2 2 2 2; 2 2 2 1 1];
    IM(:, :, 2) = [1 2 2 1 2; 2 1 2 0 2; 2 2 2 0 1; 1 1 1 0 0; 1 0 2 2 1];
    IM(:, :, 3) = [0 0 2 2 1; 0 2 1 1 2; 0 2 0 0 1; 0 2 1 0 1; 1 2 1 0 0];

    LB(:,:,1) = [-0.1 -0.2 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0];
    LB(:,:,2) = [-0.1 -0.15 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0];
    LB(:,:,3) = LB(:,:,2);

    UB(:,:,1) = [0.1 0.2 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0];
    UB(:,:,2) = [0.1 0.15 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0];
    UB(:,:,3) = UB(:,:,2);

    input_attacked = ImageStar(IM, LB, UB);
    L = AveragePooling2DLayer([3 3], [2 2], [0 0 0 0]);
    Y = L.reach_star_single_input(input_attacked);

    % ASSERTION 7: Reach with attacked image works
    assert(~isempty(Y), 'reach_star_single_input with attacked image should produce non-empty result');

    % Save regression data
    data = struct();
    data.averageMap = averageMap;
    data.y_shape = size(y);
    data.output_image_height = output_image.height;
    data.output_image_width = output_image.width;
    save_test_data(data, 'test_layers_averagePooling', 'results', 'subdir', 'nn');
end
