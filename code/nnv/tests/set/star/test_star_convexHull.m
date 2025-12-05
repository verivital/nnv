function test_star_convexHull()
    % TEST_STAR_CONVEXHULL - Test Star.convexHull() method
    %
    % Tests that:
    %   1. Convex hull produces valid output star
    %   2. Hull contains all points from both input stars

    % Create first star set
    I1 = ExamplePoly.randVrep;
    V = [0 0; 1 0; 0 1];
    I1 = Star(V', I1.A, I1.b);

    % Create second star set
    I2 = ExamplePoly.randVrep;
    V = [1 1; 1 0; 0 1];
    I2 = Star(V', I2.A, I2.b);

    % Compute convex hull
    I21 = I2.convexHull(I1);

    % ASSERTION 1: Output star is valid (non-empty)
    assert(~I21.isEmptySet, 'Convex hull should not be empty');

    % ASSERTION 2: Hull should contain points from both stars
    tol = 1e-6;
    num_samples = 15;

    % Sample from I1 and verify containment in hull bounds
    samples_I1 = I1.sample(num_samples);
    [lb_hull, ub_hull] = I21.getRanges();

    for i = 1:size(samples_I1, 2)
        s = samples_I1(:, i);
        assert(all(s >= lb_hull - tol) && all(s <= ub_hull + tol), ...
            sprintf('I1 sample %d should be in convex hull bounds', i));
    end

    % Sample from I2 and verify containment in hull bounds
    samples_I2 = I2.sample(num_samples);
    for i = 1:size(samples_I2, 2)
        s = samples_I2(:, i);
        assert(all(s >= lb_hull - tol) && all(s <= ub_hull + tol), ...
            sprintf('I2 sample %d should be in convex hull bounds', i));
    end

    % ASSERTION 3: Hull bounds should encompass both input bounds
    [lb1, ub1] = I1.getRanges();
    [lb2, ub2] = I2.getRanges();
    expected_lb = min(lb1, lb2);
    expected_ub = max(ub1, ub2);

    assert(all(lb_hull <= expected_lb + tol), ...
        'Convex hull lower bound should be <= min of input lower bounds');
    assert(all(ub_hull >= expected_ub - tol), ...
        'Convex hull upper bound should be >= max of input upper bounds');

    % Create visualization
    fig = figure;
    Star.plot(I21);
    hold on;
    Star.plot(I1);
    hold on;
    Star.plot(I2);
    legend('Convex Hull', 'I1', 'I2');
    title('Star Convex Hull Test');

    save_test_figure(fig, 'test_star_convexHull', 'convexHull', 1, 'subdir', 'set/star');

    % Save regression data
    data = struct();
    data.I1_V = I1.V;
    data.I2_V = I2.V;
    data.I21_V = I21.V;
    data.hull_lb = lb_hull;
    data.hull_ub = ub_hull;
    save_test_data(data, 'test_star_convexHull', 'results', 'subdir', 'set');
end
