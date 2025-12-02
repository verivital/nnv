function test_star_getBoxFast()
    % TEST_STAR_GETBOXFAST - Test Star.getBoxFast() vs Star.getBox()
    %
    % Tests that:
    %   1. getBoxFast() produces equivalent bounds to getBox()
    %   2. The box bounds over-approximate the star (all star samples within bounds)

    % Create a star set
    I1 = ExamplePoly.randVrep;
    V = [0 0; 1 1; 0.5 1];
    I1 = Star(V', I1.A, I1.b);

    % Compute bounding box using both methods
    B1 = I1.getBox();       % Standard method (LP-based)
    B2 = I1.getBoxFast();   % Fast method (zonotope over-approximation)

    % ASSERTION 1: Fast method should give valid box (may be larger, never smaller)
    % B2 should over-approximate B1 (or be equal)
    tol = 1e-6;
    assert(all(B2.lb <= B1.lb + tol), ...
        'getBoxFast lower bounds should be <= getBox lower bounds');
    assert(all(B2.ub >= B1.ub - tol), ...
        'getBoxFast upper bounds should be >= getBox upper bounds');

    % ASSERTION 2: Sample points from star should be within BOTH boxes
    num_samples = 20;
    samples = I1.sample(num_samples);

    for i = 1:size(samples, 2)
        s = samples(:, i);

        % Check against getBox bounds
        assert(all(s >= B1.lb - tol) && all(s <= B1.ub + tol), ...
            sprintf('Sample %d should be within getBox bounds', i));

        % Check against getBoxFast bounds (should always pass if B1 passes)
        assert(all(s >= B2.lb - tol) && all(s <= B2.ub + tol), ...
            sprintf('Sample %d should be within getBoxFast bounds', i));
    end

    % Create visualization
    fig = figure;
    Box.plot(B2);
    hold on;
    Box.plot(B1);
    hold on;
    Star.plot(I1);
    legend('getBoxFast (B2)', 'getBox (B1)', 'Star I1');
    title('Star Bounding Box Comparison');
    xlabel('x_1'); ylabel('x_2');

    % Save figure
    save_test_figure(fig, 'test_star_getBoxFast', 'boxComparison', 1, 'subdir', 'set/star');

    % Save regression data
    data = struct();
    data.B1_lb = B1.lb;
    data.B1_ub = B1.ub;
    data.B2_lb = B2.lb;
    data.B2_ub = B2.ub;
    data.I1_V = I1.V;
    data.I1_C = I1.C;
    data.I1_d = I1.d;
    save_test_data(data, 'test_star_getBoxFast', 'results', 'subdir', 'set');
end
