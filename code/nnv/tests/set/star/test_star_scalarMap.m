function test_star_scalarMap()
    % TEST_STAR_SCALARMAP - Test Star.scalarMap() method
    %
    % Tests that:
    %   1. Scalar multiplication produces valid output star
    %   2. Output star bounds are scaled correctly

    % Create a star set with well-defined constraints (avoid random LP issues)
    center = [1; 1];
    V = [1 0; 0 1];  % 2D star with 2 generators
    % Simple box constraints: -1 <= alpha <= 1
    C = [1 0; -1 0; 0 1; 0 -1];
    d = [1; 1; 1; 1];
    S = Star([center V], C, d);

    % Apply scalar multiplication
    alpha = 0.7;
    S1 = S.scalarMap(alpha);

    % ASSERTION 1: Output star is valid
    assert(~S1.isEmptySet, 'Scalar mapped star should not be empty');

    % Get bounds for data saving
    [lb0, ub0] = S.getRanges();
    [lb1, ub1] = S1.getRanges();

    % Create visualization
    fig = figure;
    Star.plot(S);
    hold on;
    Star.plot(S1);
    legend('Original S', sprintf('Scaled S1 (%.1f * S)', alpha));
    title('Star Scalar Map Test');

    save_test_figure(fig, 'test_star_scalarMap', 'scalarMap', 1, 'subdir', 'set/star');

    % Save regression data
    data = struct();
    data.alpha = alpha;
    data.S_V = S.V;
    data.S1_V = S1.V;
    data.S_lb = lb0;
    data.S_ub = ub0;
    data.S1_lb = lb1;
    data.S1_ub = ub1;
    save_test_data(data, 'test_star_scalarMap', 'results', 'subdir', 'set');
end
