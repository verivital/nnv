function test_star_getBox()
    % TEST_STAR_GETBOX - Test Star.getBox() method
    %
    % Tests that:
    %   1. getBox() returns correct bounds for a box-derived star
    %   2. Returned box contains all points in the star

    % Create a box and convert to star
    lb = [49; 25; 9; 20];
    ub = [51; 25.2; 11; 20.2];
    B1 = Box(lb, ub);
    init_set = B1.toStar();

    % Get box from star
    tic;
    B2 = init_set.getBox();
    elapsed = toc;

    % ASSERTION 1: Recovered bounds should match original box
    tol = 1e-6;
    assert(all(abs(B2.lb - lb) < tol), 'Lower bounds should match original box');
    assert(all(abs(B2.ub - ub) < tol), 'Upper bounds should match original box');

    % ASSERTION 2: All sampled points should be within box
    num_samples = 20;
    samples = init_set.sample(num_samples);
    for i = 1:size(samples, 2)
        s = samples(:, i);
        assert(all(s >= B2.lb - tol) && all(s <= B2.ub + tol), ...
            sprintf('Sample %d should be within box bounds', i));
    end

    % ASSERTION 3: Performance check (should be fast for simple case)
    assert(elapsed < 5.0, 'getBox should complete in under 5 seconds');

    % Save regression data
    data = struct();
    data.original_lb = lb;
    data.original_ub = ub;
    data.recovered_lb = B2.lb;
    data.recovered_ub = B2.ub;
    data.elapsed_time = elapsed;
    save_test_data(data, 'test_star_getBox', 'results', 'subdir', 'set');
end
