function test_star_MinkowskiSum()
    % TEST_STAR_MINKOWSKISUM - Test Star.MinkowskiSum() method
    %
    % Tests that:
    %   1. Minkowski sum produces valid output star
    %   2. Sum is sound: for all x in S1, y in S2, x+y is in S1+S2

    % Create first star set
    V = [1 1 0; 0 1 0; 0 0 1];
    C = [1 0; -1 0; 0 1; 0 -1];
    d = [1; 1; 1; 1];  % -1 <= a[1] <= 1, -1 <= a[2] <= 1
    S1 = Star(V, C, d);

    % Create second star set
    V2 = [1 0; 0 1; 1 1];
    C2 = [1; -1];
    d2 = [0.5; 0.5];
    S2 = Star(V2, C2, d2);

    % Create third star via affine map
    W = [2 1 1; 1 0 2; 0 1 0];
    b = [0.5; 0.5; 0];
    S3 = S1.affineMap(W, b);

    % Compute Minkowski sums
    S12 = S1.MinkowskiSum(S2);
    S13 = S1.MinkowskiSum(S3);

    % ASSERTION 1: Output stars are valid (non-empty)
    assert(~S12.isEmptySet, 'Minkowski sum S12 should not be empty');
    assert(~S13.isEmptySet, 'Minkowski sum S13 should not be empty');

    % ASSERTION 2: Verify soundness - sampled sums should be in result
    num_samples = 10;
    tol = 1e-6;

    % Sample from S1 and S2
    samples_S1 = S1.sample(num_samples);
    samples_S2 = S2.sample(num_samples);

    for i = 1:min(size(samples_S1, 2), size(samples_S2, 2))
        x = samples_S1(:, i);
        y = samples_S2(:, i);
        sum_xy = x + y;

        % Check containment in S12 bounds
        [lb, ub] = S12.getRanges();
        assert(all(sum_xy >= lb - tol) && all(sum_xy <= ub + tol), ...
            sprintf('Sum point %d should be in Minkowski sum bounds', i));
    end

    % Create visualization - figure 1
    fig1 = figure;
    Star.plot(S12);
    hold on;
    Star.plot(S1);
    hold on;
    Star.plot(S2);
    legend('S1 + S2', 'S1', 'S2');
    title('Minkowski Sum: S12 = S1 + S2');

    save_test_figure(fig1, 'test_star_MinkowskiSum', 'sum_S12', 1, 'subdir', 'set/star');

    % Create visualization - figure 2
    fig2 = figure;
    Star.plot(S13);
    title('Minkowski Sum: S13 = S1 + S3');

    save_test_figure(fig2, 'test_star_MinkowskiSum', 'sum_S13', 2, 'subdir', 'set/star');

    % Save regression data
    data = struct();
    data.S1_V = S1.V;
    data.S2_V = S2.V;
    data.S12_V = S12.V;
    [data.S12_lb, data.S12_ub] = S12.getRanges();
    [data.S13_lb, data.S13_ub] = S13.getRanges();
    save_test_data(data, 'test_star_MinkowskiSum', 'results', 'subdir', 'set');
end
