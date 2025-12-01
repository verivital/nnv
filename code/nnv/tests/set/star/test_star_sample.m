function test_star_sample()
    % TEST_STAR_SAMPLE - Test Star.sample() method
    %
    % Tests that:
    %   1. sample() returns points inside the star
    %   2. All sampled points satisfy the star constraints

    % Create a star set
    I = ExamplePoly.randVrep;
    V = [0 0; 1 0; 0 1];
    S = Star(V', I.A, I.b);

    % Sample points from star
    num_samples = 100;
    samples = S.sample(num_samples);

    % ASSERTION 1: Some samples should be returned
    assert(~isempty(samples), 'sample() should return non-empty result');
    assert(size(samples, 2) > 0, 'At least one sample should be valid');

    % ASSERTION 2: All samples should be within star bounds
    tol = 1e-6;
    [lb, ub] = S.getRanges();
    for i = 1:size(samples, 2)
        s = samples(:, i);
        assert(all(s >= lb - tol) && all(s <= ub + tol), ...
            sprintf('Sample %d should be within star bounds', i));
    end

    % ASSERTION 3: All samples should satisfy star constraints
    % For a point to be in the star: x = c + V*alpha where C*alpha <= d
    for i = 1:size(samples, 2)
        s = samples(:, i);
        % Check containment
        contained = S.contains(s);
        assert(contained == 1, sprintf('Sample %d should be contained in star', i));
    end

    % Create visualization
    fig = figure;
    Star.plot(S);
    hold on;
    plot(samples(1, :), samples(2, :), 'r*', 'MarkerSize', 4);
    legend('Star S', 'Samples');
    title(sprintf('Star Sampling Test (%d samples)', size(samples, 2)));

    save_test_figure(fig, 'test_star_sample', 'samples', 1, 'subdir', 'set/star');

    % Save regression data
    data = struct();
    data.num_requested = num_samples;
    data.num_returned = size(samples, 2);
    data.star_lb = lb;
    data.star_ub = ub;
    data.S_V = S.V;
    save_test_data(data, 'test_star_sample', 'results', 'subdir', 'set');
end
