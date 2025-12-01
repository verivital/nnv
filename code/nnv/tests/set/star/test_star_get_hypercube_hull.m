function test_star_get_hypercube_hull()
    % TEST_STAR_GET_HYPERCUBE_HULL - Test Star.get_hypercube_hull() static method
    %
    % Tests that:
    %   1. Hypercube (box) hull of array of stars is computed correctly
    %   2. Box contains all points from input stars

    % Create first star set
    I1 = ExamplePoly.randVrep;
    V = [0 0; 1 0; 0 1];
    I1 = Star(V', I1.A, I1.b);

    % Create second star set
    I2 = ExamplePoly.randVrep;
    V = [1 1; 1 0; 0 1];
    I2 = Star(V', I2.A, I2.b);

    % Compute hypercube hull
    I3 = Star.get_hypercube_hull([I1 I2]);

    % ASSERTION 1: Result is a valid box
    assert(~isempty(I3), 'get_hypercube_hull should return a valid box');
    assert(isa(I3, 'Box'), 'Result should be a Box object');

    % ASSERTION 2: Box should contain all points from input stars
    tol = 1e-6;
    num_samples = 15;

    % Check I1 samples
    samples_I1 = I1.sample(num_samples);
    for i = 1:size(samples_I1, 2)
        s = samples_I1(:, i);
        assert(all(s >= I3.lb - tol) && all(s <= I3.ub + tol), ...
            sprintf('I1 sample %d should be in hypercube hull', i));
    end

    % Check I2 samples
    samples_I2 = I2.sample(num_samples);
    for i = 1:size(samples_I2, 2)
        s = samples_I2(:, i);
        assert(all(s >= I3.lb - tol) && all(s <= I3.ub + tol), ...
            sprintf('I2 sample %d should be in hypercube hull', i));
    end

    % ASSERTION 3: Box bounds should encompass both input bounds
    [lb1, ub1] = I1.getRanges();
    [lb2, ub2] = I2.getRanges();
    expected_lb = min(lb1, lb2);
    expected_ub = max(ub1, ub2);

    assert(all(I3.lb <= expected_lb + tol), ...
        'Box lower bounds should be <= min of input lower bounds');
    assert(all(I3.ub >= expected_ub - tol), ...
        'Box upper bounds should be >= max of input upper bounds');

    % Create visualization
    fig = figure;
    Box.plot(I3);
    hold on;
    Star.plot(I2);
    hold on;
    Star.plot(I1);
    legend('Hypercube Hull (Box)', 'I2', 'I1');
    title('Star.get\_hypercube\_hull() Test');

    save_test_figure(fig, 'test_star_get_hypercube_hull', 'hull', 1, 'subdir', 'set/star');

    % Save regression data
    data = struct();
    data.I1_V = I1.V;
    data.I2_V = I2.V;
    data.I3_lb = I3.lb;
    data.I3_ub = I3.ub;
    save_test_data(data, 'test_star_get_hypercube_hull', 'results', 'subdir', 'set');
end
