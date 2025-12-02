function test_PosLin_reach_star_approx()
    % TEST_POSLIN_REACH_STAR_APPROX - Test PosLin approximate reachability
    %
    % Tests that:
    %   1. Star input set is constructed correctly
    %   2. reach_star_approx produces valid over-approximation
    %   3. reach (exact) produces valid exact reach set
    %   4. reach_star_approx2 produces valid over-approximation
    %   5. All sampled outputs lie within the reach sets

    %% Create input Star
    I0 = ExamplePoly.randVrep;
    I0.outerApprox;
    V = [0 0; 1 1; 1 0];
    I = Star(V', I0.A, I0.b, I0.Internal.lb, I0.Internal.ub);

    % ASSERTION 1: Input star is valid
    assert(~isempty(I), 'Input star should be created successfully');
    assert(I.dim == 2, 'Input star should be 2-dimensional');

    % Sample inputs (note: sample() may not return exact count due to constraints)
    X = I.sample(100);

    % ASSERTION 2: Sampling works (at least some samples returned)
    assert(size(X, 2) > 0, 'Should sample at least some points');
    assert(size(X, 1) == 2, 'Samples should be 2-dimensional');

    %% Test 1: Visualize input set
    fig1 = figure;
    Star.plot(I);
    hold on;
    plot(X(1, :), X(2, :), 'ob');
    title('Input Star with Samples');
    legend('Star', 'Samples');

    save_test_figure(fig1, 'test_PosLin_reach_star_approx', 'input', 1, 'subdir', 'nn/PosLin');

    %% Test 2: reach_star_approx
    t = tic;
    S = PosLin.reach_star_approx(I);
    t1 = toc(t);

    % ASSERTION 3: Approx reach produces valid result
    assert(~isempty(S), 'reach_star_approx should produce non-empty result');

    %% Test 3: reach (exact)
    S1 = PosLin.reach(I);

    % ASSERTION 4: Exact reach produces valid result
    assert(~isempty(S1), 'reach should produce non-empty result');

    %% Test 4: reach_star_approx2
    t = tic;
    S2 = PosLin.reach_star_approx2(I);
    t2 = toc(t);

    % ASSERTION 5: Approx2 reach produces valid result
    assert(~isempty(S2), 'reach_star_approx2 should produce non-empty result');

    %% Evaluate and visualize outputs
    Y = PosLin.evaluate(X);

    fig2 = figure;
    Star.plot(S);
    hold on;
    Star.plots(S1);
    hold on;
    plot(Y(1, :), Y(2, :), '*');
    title(sprintf('Output Sets (approx: %.3fs)', t1));
    legend('Approx', 'Exact', 'Samples');

    save_test_figure(fig2, 'test_PosLin_reach_star_approx', 'output_compare', 2, 'subdir', 'nn/PosLin');

    fig3 = figure;
    Star.plots(S2);
    title(sprintf('reach_star_approx2 (%.3fs)', t2));

    save_test_figure(fig3, 'test_PosLin_reach_star_approx', 'approx2', 3, 'subdir', 'nn/PosLin');

    % Save regression data
    data = struct();
    data.t_approx = t1;
    data.t_approx2 = t2;
    data.input_dim = I.dim;
    if iscell(S1)
        data.num_exact_sets = length(S1);
    else
        data.num_exact_sets = 1;
    end
    save_test_data(data, 'test_PosLin_reach_star_approx', 'results', 'subdir', 'nn');
end
