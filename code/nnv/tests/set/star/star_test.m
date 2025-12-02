function star_test()
    % STAR_TEST - Basic Star class functionality test
    %
    % Tests that:
    %   1. Star construction works correctly
    %   2. Affine map transforms stars correctly
    %   3. Minkowski sum is computed correctly
    %   4. getBox returns valid bounding box

    % Create first star set
    V = [1 1 0; 0 1 0; 0 0 1];
    C = [1 0; -1 0; 0 1; 0 -1];
    d = [1; 1; 1; 1];  % -1 <= a[1] <= 1, -1 <= a[2] <= 1
    S1 = Star(V, C, d);

    % ASSERTION 1: Star is valid
    assert(~S1.isEmptySet, 'Star S1 should not be empty');
    assert(S1.dim == 3, 'Star S1 should have dimension 3');

    % Apply affine transformation
    W = [2 1 1; 1 0 2];
    b = [0.5; 0.5];
    S11 = S1.affineMap(W, b);

    % ASSERTION 2: Affine map result is valid
    assert(~S11.isEmptySet, 'Transformed star S11 should not be empty');
    assert(S11.dim == 2, 'Transformed star S11 should have dimension 2');

    % Create second star set
    V2 = [1 0; 0 1; 1 1];
    C2 = [1; -1];
    d2 = [0.5; 0.5];
    S2 = Star(V2, C2, d2);

    % Compute Minkowski sum
    S3 = S1.MinkowskiSum(S2);

    % ASSERTION 3: Minkowski sum is valid and contains sum of samples
    assert(~S3.isEmptySet, 'Minkowski sum S3 should not be empty');

    tol = 1e-6;
    samples_S1 = S1.sample(5);
    samples_S2 = S2.sample(5);
    [lb3, ub3] = S3.getRanges();

    for i = 1:min(size(samples_S1, 2), size(samples_S2, 2))
        sum_point = samples_S1(:, i) + samples_S2(:, i);
        assert(all(sum_point >= lb3 - tol) && all(sum_point <= ub3 + tol), ...
            sprintf('Sum point %d should be in Minkowski sum bounds', i));
    end

    % Get bounding box of S2
    B2 = S2.getBox();

    % ASSERTION 4: Box bounds contain star samples
    samples_S2_check = S2.sample(10);
    for i = 1:size(samples_S2_check, 2)
        s = samples_S2_check(:, i);
        assert(all(s >= B2.lb - tol) && all(s <= B2.ub + tol), ...
            sprintf('S2 sample %d should be within box bounds', i));
    end

    % Create visualization
    fig = figure;
    Star.plot(S1);
    hold on;
    Star.plot(S2);
    hold on;
    Star.plot(S3);
    legend('S1', 'S2', 'S3 = S1 + S2');
    title('Star Test: Construction, Transform, Minkowski Sum');

    save_test_figure(fig, 'star_test', 'basic', 1, 'subdir', 'set/star');

    % Save regression data
    data = struct();
    data.S1_V = S1.V;
    data.S2_V = S2.V;
    data.S3_V = S3.V;
    data.S11_dim = S11.dim;
    data.B2_lb = B2.lb;
    data.B2_ub = B2.ub;
    save_test_data(data, 'star_test', 'results', 'subdir', 'set');
end
