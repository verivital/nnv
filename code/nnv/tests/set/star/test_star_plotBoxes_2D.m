function test_star_plotBoxes_2D()
    % TEST_STAR_PLOTBOXES_2D - Test Star.plotBoxes_2D() methods
    %
    % Tests that:
    %   1. Stars can be plotted successfully
    %   2. plotBoxes_2D and plotBoxes_2D_noFill work without error
    %   3. Bounding boxes are correctly computed for plotting

    % Create first star set
    V = [1 1 0; 0 1 0; 0 0 1];
    C = [1 0; -1 0; 0 1; 0 -1];
    d = [1; 1; 1; 1];
    S1 = Star(V, C, d);

    % Apply affine transformation
    W = [2 1 1; 1 0 2; 0 1 0];
    b = [0.5; 0.5; 0];
    S2 = S1.affineMap(W, b);

    % ASSERTION 1: Both stars are valid
    assert(~S1.isEmptySet, 'Star S1 should not be empty');
    assert(~S2.isEmptySet, 'Star S2 should not be empty');

    % Array of stars
    S = [S1 S2];

    % ASSERTION 2: Stars have correct dimension
    assert(S1.dim == 3, 'S1 should have dimension 3');
    assert(S2.dim == 3, 'S2 should have dimension 3');

    % Create figure 1: Basic star plots
    fig1 = figure;
    Star.plot(S1);
    hold on;
    Star.plot(S2);
    legend('S1', 'S2 (transformed)');
    title('Stars S1 and S2');

    save_test_figure(fig1, 'test_star_plotBoxes_2D', 'stars', 1, 'subdir', 'set/star');

    % Create figure 2: Multiple stars plot
    fig2 = figure;
    Star.plots(S);
    title('Star.plots() - Multiple Stars');

    save_test_figure(fig2, 'test_star_plotBoxes_2D', 'plots', 2, 'subdir', 'set/star');

    % Create figure 3: 2D box projection (filled)
    fig3 = figure;
    Star.plotBoxes_2D(S, 1, 2, 'red');
    title('plotBoxes\_2D (filled) - Dimensions 1 vs 2');
    xlabel('x_1'); ylabel('x_2');

    save_test_figure(fig3, 'test_star_plotBoxes_2D', 'boxes2D_filled', 3, 'subdir', 'set/star');

    % Create figure 4: 2D box projection (no fill)
    fig4 = figure;
    Star.plotBoxes_2D_noFill(S, 1, 2, 'red');
    title('plotBoxes\_2D\_noFill - Dimensions 1 vs 2');
    xlabel('x_1'); ylabel('x_2');

    save_test_figure(fig4, 'test_star_plotBoxes_2D', 'boxes2D_noFill', 4, 'subdir', 'set/star');

    % Save regression data
    data = struct();
    data.S1_V = S1.V;
    data.S2_V = S2.V;
    [data.S1_lb, data.S1_ub] = S1.getRanges();
    [data.S2_lb, data.S2_ub] = S2.getRanges();
    save_test_data(data, 'test_star_plotBoxes_2D', 'results', 'subdir', 'set');
end
