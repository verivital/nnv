function test_star_plotBoxes_3D()
    % TEST_STAR_PLOTBOXES_3D - Test 3D box plotting for Star sets
    %
    % Tests that:
    %   1. Star from random H-rep polytope is constructed correctly
    %   2. 3D box plotting works without errors

    %% Test 1: Create 3D Star and plot boxes
    rng(1);  % For reproducibility
    center = [1; 1; 1];
    V = [1 0 1; 0 1 1; 1 0 0];
    P = ExamplePoly.randHrep('d', 3);
    S = Star([center V], P.A, P.b);

    % ASSERTION 1: Star is valid
    assert(~isempty(S), 'Star should be created successfully');
    assert(S.dim == 3, 'Star should be 3-dimensional');

    % ASSERTION 2: Star has correct basis
    expected_V = [center V];
    assert(size(S.V, 1) == 3, 'Star basis should have 3 rows');
    assert(size(S.V, 2) == 4, 'Star basis should have 4 columns (center + 3 generators)');

    % Plot 3D boxes
    fig1 = figure;
    Star.plotBoxes_3D(S, 1, 2, 3, 'green');
    title('Star 3D Boxes');

    % ASSERTION 3: Figure was created successfully
    assert(ishandle(fig1), 'Figure should be created successfully');

    save_test_figure(fig1, 'test_star_plotBoxes_3D', 'boxes_3d', 1, 'subdir', 'set/star');

    % Save regression data
    data = struct();
    data.center = center;
    data.V = V;
    data.star_dim = S.dim;
    data.star_nVar = S.nVar;
    save_test_data(data, 'test_star_plotBoxes_3D', 'results', 'subdir', 'set');
end
