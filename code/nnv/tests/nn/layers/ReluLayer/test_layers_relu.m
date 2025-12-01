function test_layers_relu()
    % TEST_LAYERS_RELU - Test ReluLayer operations
    %
    % Tests that:
    %   1. ReluLayer constructor works correctly
    %   2. ReluLayer evaluate produces correct output
    %   3. ReluLayer reach (approx-star) works correctly
    %   4. ReluLayer reach (exact-star) works correctly
    %   5. ReluLayer reach (zono) works correctly

    %% Test 1: ReluLayer constructor
    L = ReluLayer();
    L1 = ReluLayer('relu1');

    % ASSERTION 1: Constructors work
    assert(~isempty(L), 'ReluLayer default constructor should work');
    assert(~isempty(L1), 'ReluLayer named constructor should work');

    %% Test 2: ReluLayer evaluate
    rl = ReluLayer();
    % Image input set
    IM(:,:,1) = [-1 1 0 -1; 0 0 1 -1; 1 0 -1 0; 1 -1 -1 1];
    IM(:,:,2) = [0 1 0 0; 1 0 0 -1; 0 1 -1 0; 0 1 0 -1];
    IM(:,:,3) = [1 -1 1 1; 1 -1 0 1; 0 1 -1 0; 1 0 -1 0];

    output = rl.evaluate(IM);

    IM_out(:,:,1) = [0 1 0 0; 0 0 1 0; 1 0 0 0; 1 0 0 1];
    IM_out(:,:,2) = [0 1 0 0; 1 0 0 0; 0 1 0 0; 0 1 0 0];
    IM_out(:,:,3) = [1 0 1 1; 1 0 0 1; 0 1 0 0; 1 0 0 0];

    % ASSERTION 2: Evaluate produces correct output
    assert(isequal(output, IM_out), 'ReluLayer evaluate should produce correct output');

    %% Test 3: ReluLayer reach star approx
    L = ReluLayer();
    IM(:,:,1) = [1 1 0 1; 0 0 1 1; 1 0 1 0; 1 1 1 1];
    IM(:,:,2) = [0 1 0 0; 1 0 0 1; 0 1 1 0; 0 0 0 1];
    IM(:,:,3) = [1 1 1 1; 1 1 0 1; 0 1 1 0; 1 0 1 0];

    LB(:,:,1) = [-0.1 -0.2 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
    LB(:,:,2) = [-0.1 -0.15 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
    LB(:,:,3) = LB(:,:,2);

    UB(:,:,1) = [0.1 0.2 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
    UB(:,:,2) = [0.1 0.15 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
    UB(:,:,3) = UB(:,:,2);

    image = ImageStar(IM, LB, UB);

    tic;
    images_approx = L.reach(image, 'approx-star');
    t_approx = toc;

    % ASSERTION 3: Approx reach produces valid output
    assert(~isempty(images_approx), 'ReluLayer approx-star reach should produce non-empty result');

    %% Test 4: ReluLayer reach star exact
    L = ReluLayer();
    V(:, :, 1, 1) = [-1 1; 0 2];
    basis = zeros(2, 2);
    basis(1, 1) = 1;
    basis(2, 1) = 1;
    V(:, :, 1, 2) = basis;

    C = [1; -1];
    d = [2; 2];
    pred_lb = -2;
    pred_ub = 2;

    in_image = ImageStar(V, C, d, pred_lb, pred_ub);

    images_exact = L.reach(in_image, 'exact-star');

    % ASSERTION 4: Exact reach produces valid output
    assert(~isempty(images_exact), 'ReluLayer exact-star reach should produce non-empty result');

    %% Test 5: ReluLayer reach zono
    L = ReluLayer();
    LB(:,:,1) = [-0.1 -0.2 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
    LB(:,:,2) = [-0.1 -0.15 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
    LB(:,:,3) = LB(:,:,2);

    UB(:,:,1) = [0 0.2 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
    UB(:,:,2) = [0.1 0.15 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
    UB(:,:,3) = UB(:,:,2);

    image_zono = ImageZono(LB, UB);

    image_out = L.reach_zono(image_zono);

    % ASSERTION 5: Zono reach produces valid output
    assert(~isempty(image_out), 'ReluLayer zono reach should produce non-empty result');

    % Save regression data
    data = struct();
    data.t_approx = t_approx;
    data.IM = IM;
    data.IM_out = output;
    if iscell(images_approx)
        data.num_approx_images = length(images_approx);
    else
        data.num_approx_images = 1;
    end
    if iscell(images_exact)
        data.num_exact_images = length(images_exact);
    else
        data.num_exact_images = 1;
    end
    save_test_data(data, 'test_layers_relu', 'results', 'subdir', 'nn');
end
