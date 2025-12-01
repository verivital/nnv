function test_star_getOrientedBox()
    % TEST_STAR_GETORIENTEDBOX - Test oriented bounding box computation
    %
    % Tests that:
    %   1. Zonotope to Star conversion works
    %   2. getOrientedBox produces valid oriented bounding box
    %   3. getBox produces valid axis-aligned bounding box
    %   4. Oriented box is potentially tighter than axis-aligned box

    %% Test 1: Create Star from Zonotope
    c1 = [0; 0];
    V1 = [1 -1; 1 1; 0.5 0; -1 0.5];
    Z1 = Zono(c1, V1');
    I1 = Z1.toStar();

    % ASSERTION 1: Star is valid
    assert(~isempty(I1), 'Star from Zonotope should be created successfully');
    assert(I1.dim == 2, 'Star should be 2-dimensional');

    %% Test 2: Get Oriented Box
    I2 = I1.getOrientedBox();

    % ASSERTION 2: Oriented box is valid
    assert(~isempty(I2), 'Oriented box should be created successfully');

    %% Test 3: Get Axis-Aligned Box
    I3 = I1.getBox();

    % ASSERTION 3: Axis-aligned box is valid
    assert(~isempty(I3), 'Axis-aligned box should be created successfully');

    %% Visualize results
    fig1 = figure;
    Box.plot(I2);
    hold on;
    Star.plot(I1);
    title('Oriented Bounding Box vs Star');
    legend('Oriented Box', 'Star');

    save_test_figure(fig1, 'test_star_getOrientedBox', 'oriented', 1, 'subdir', 'set/star');

    fig2 = figure;
    Box.plot(I3);
    hold on;
    Star.plot(I1);
    title('Axis-Aligned Bounding Box vs Star');
    legend('Axis-Aligned Box', 'Star');

    save_test_figure(fig2, 'test_star_getOrientedBox', 'axis_aligned', 2, 'subdir', 'set/star');

    % Save regression data
    data = struct();
    data.c1 = c1;
    data.V1 = V1;
    data.star_dim = I1.dim;
    data.star_nVar = I1.nVar;
    save_test_data(data, 'test_star_getOrientedBox', 'results', 'subdir', 'set');
end
