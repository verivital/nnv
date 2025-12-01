function test_star_intersect_halfSpace()
    % TEST_STAR_INTERSECT_HALFSPACE - Test Star.intersectHalfSpace() method
    %
    % Tests that:
    %   1. Intersection produces valid output star
    %   2. All points in result satisfy the halfspace constraint
    %   3. Result is a subset of the original star

    % Create a star set
    center = [1; 1];
    V = [1 0 1; 0 1 1];
    P = ExamplePoly.randHrep('d', 3);
    S = Star([center V], P.A, P.b);

    % Define halfspace: x[1] >= 1  =>  -x[1] <= -1
    H = [-1 0];
    g = [-1];

    % Compute intersection
    S1 = S.intersectHalfSpace(H, g);

    % ASSERTION 1: Check if intersection is empty (valid result either way)
    if S1.isEmptySet
        warning('Intersection is empty - halfspace does not overlap with star');
    else
        % ASSERTION 2: All sampled points should satisfy the halfspace constraint
        num_samples = 20;
        samples = S1.sample(num_samples);
        tol = 1e-6;

        for i = 1:size(samples, 2)
            s = samples(:, i);
            % Check H*s <= g (which is -s(1) <= -1, i.e., s(1) >= 1)
            constraint_value = H * s;
            assert(constraint_value <= g + tol, ...
                sprintf('Sample %d should satisfy halfspace constraint (H*x <= g)', i));
        end

        % ASSERTION 3: Result bounds should be within original bounds
        [lb_orig, ub_orig] = S.getRanges();
        [lb_int, ub_int] = S1.getRanges();

        assert(all(lb_int >= lb_orig - tol), ...
            'Intersection lower bounds should be >= original lower bounds');
        assert(all(ub_int <= ub_orig + tol), ...
            'Intersection upper bounds should be <= original upper bounds');

        % ASSERTION 4: First dimension should be >= 1 after intersection
        assert(lb_int(1) >= 1 - tol, ...
            'After intersection with x1 >= 1, lower bound of x1 should be >= 1');
    end

    % Create visualization
    fig = figure;
    Star.plot(S);
    hold on;
    if ~S1.isEmptySet
        Star.plot(S1);
        legend('Original S', 'Intersection S1 (x_1 >= 1)');
    else
        legend('Original S (intersection empty)');
    end
    title('Star Halfspace Intersection Test');

    save_test_figure(fig, 'test_star_intersect_halfSpace', 'intersection', 1, 'subdir', 'set/star');

    % Save regression data
    data = struct();
    data.H = H;
    data.g = g;
    data.S_V = S.V;
    data.S1_V = S1.V;
    data.S1_C = S1.C;
    data.S1_d = S1.d;
    data.S1_isEmpty = S1.isEmptySet;
    if ~S1.isEmptySet
        [data.S1_lb, data.S1_ub] = S1.getRanges();
    end
    save_test_data(data, 'test_star_intersect_halfSpace', 'results', 'subdir', 'set');
end
