function test_star_affineMap()
    % TEST_STAR_AFFINEMAP - Test Star.affineMap() method
    %
    % Tests that:
    %   1. Affine transformation produces valid output star
    %   2. Transformation is sound (all input points map to output set)

    % Create a star set
    center = [1; 1];
    V = [1 0 1; 0 1 1];
    P = ExamplePoly.randHrep('d', 3);
    V = [center V];
    C = P.A;
    d = P.b;
    S = Star(V, C, d);

    % Define affine transformation: y = W*x + b
    W = [1 -1; 1 1];
    b = [0.5; 0.5];

    % Apply transformation
    S1 = S.affineMap(W, b);

    % ASSERTION 1: Output star is valid (non-empty)
    assert(~S1.isEmptySet, 'Transformed star S1 should not be empty');

    % ASSERTION 2: Verify transformation is sound using sampling
    % Sample points from input star and verify they map into output star
    num_samples = 20;
    samples = S.sample(num_samples);

    for i = 1:size(samples, 2)
        x = samples(:, i);
        y = W * x + b;  % Expected output

        % Check containment in output star
        contained = S1.contains(y);
        assert(contained == 1, sprintf('Transformed point %d should be in output star', i));
    end

    % Create visualization
    fig = figure;
    subplot(1, 2, 1);
    Star.plot(S);
    title('Input Star S');
    xlabel('x_1'); ylabel('x_2');

    subplot(1, 2, 2);
    Star.plot(S1);
    title('Transformed Star S1 = W*S + b');
    xlabel('y_1'); ylabel('y_2');

    sgtitle('Star Affine Map Test');

    % Save figure
    save_test_figure(fig, 'test_star_affineMap', 'affineMap', 1, 'subdir', 'set/star');

    % Test with Box-defined star
    lb = [0; -1; 0];
    ub = [1; 1; 1];
    S3 = Star(lb, ub);

    % Use fixed W for reproducibility in regression tests
    rng(42);
    W2 = rand(3, 3);
    S4 = S3.affineMap(W2, []);

    % ASSERTION 3: Box-to-Star transformation is valid
    assert(~S4.isEmptySet, 'Box-derived star transformation should produce non-empty set');

    % Verify soundness for box transformation
    for i = 1:10
        x = lb + (ub - lb) .* rand(3, 1);
        y = W2 * x;
        [lb_out, ub_out] = S4.getRanges();
        assert(all(y >= lb_out - 1e-6) && all(y <= ub_out + 1e-6), ...
            sprintf('Transformed box point %d should be in output bounds', i));
    end

    % Save regression data
    data = struct();
    data.W = W;
    data.b = b;
    data.S_V = S.V;
    data.S_C = S.C;
    data.S_d = S.d;
    data.S1_V = S1.V;
    data.S1_C = S1.C;
    data.S1_d = S1.d;
    data.W2 = W2;
    [data.S4_lb, data.S4_ub] = S4.getRanges();
    save_test_data(data, 'test_star_affineMap', 'results', 'subdir', 'set');
end
