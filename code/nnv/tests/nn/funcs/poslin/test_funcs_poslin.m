function test_funcs_poslin()
    % TEST_FUNCS_POSLIN - Test PosLin (ReLU) reachability methods
    %
    % Tests that:
    %   1. PosLin.evaluate computes correct output
    %   2. Exact and approximate star methods work correctly
    %   3. Zonotope methods work correctly
    %   4. Abstract domain methods work correctly

    % Create random set
    I = ExamplePoly.randVrep;
    I.outerApprox;
    V = [0 0; 1 0; 0 1];
    I_poly = Star(V', I.A, I.b, I.Internal.lb, I.Internal.ub);

    lb = [-0.5; -0.5];
    ub = [0.5; 0.5];
    B = Box(lb, ub);
    I_zono = B.toZono;
    A = [0.5 1; 1.5 -2];
    I_zono = I_zono.affineMap(A, []);
    I_star = I_zono.toStar;

    %% Test 1: PosLin evaluate
    x = [-1; 0.5; 2];
    y = PosLin.evaluate(x);

    % ASSERTION 1: Evaluate produces correct output
    assert(isequal(y, [0; 0.5; 2]), 'PosLin.evaluate should return [0; 0.5; 2]');

    %% Test 2: PosLin reach exact star
    S_exact = PosLin.reach(I_star, 'exact-star');

    % ASSERTION 2: Exact reach produces valid output
    assert(~isempty(S_exact), 'Exact-star reach should produce non-empty result');

    %% Test 3: PosLin reach approx star
    S_approx = PosLin.reach(I_star, 'approx-star');

    % ASSERTION 3: Approx reach produces valid output
    assert(~isempty(S_approx), 'Approx-star reach should produce non-empty result');

    %% Test 4: PosLin reach approx zono
    Z = PosLin.reach(I_zono, 'approx-zono');

    % ASSERTION 4: Zono reach produces valid output
    assert(~isempty(Z), 'Approx-zono reach should produce non-empty result');
    assert(isa(Z, 'Zono'), 'Approx-zono result should be a Zono');

    %% Test 5: PosLin reach abstract domain
    S_abs = PosLin.reach_abstract_domain(I_poly);

    % ASSERTION 5: Abstract domain reach produces valid output
    assert(~isempty(S_abs), 'Abstract domain reach should produce non-empty result');

    %% Test 6: PosLin reach star approx
    S_approx2 = PosLin.reach_star_approx(I_poly);

    % ASSERTION 6: Star approx reach produces valid output
    assert(~isempty(S_approx2), 'Star approx reach should produce non-empty result');

    %% Test 7: PosLin reach default (exact)
    S_default = PosLin.reach(I_poly);

    % ASSERTION 7: Default reach produces valid output
    assert(~isempty(S_default), 'Default reach should produce non-empty result');

    %% Test 8: PosLin reach zono approx
    Z_approx = PosLin.reach_zono_approx(I_zono);

    % ASSERTION 8: Zono approx reach produces valid output
    assert(~isempty(Z_approx), 'Zono approx reach should produce non-empty result');

    %% Test 9: PosLin stepReach Star Approx
    S_step = PosLin.stepReachStarApprox(I_star, 1);

    % ASSERTION 9: Step reach produces valid output
    assert(~isempty(S_step), 'Step reach should produce non-empty result');

    % Create visualization
    fig = figure;
    subplot(2, 2, 1);
    Star.plot(I_star);
    title('Input Star');

    subplot(2, 2, 2);
    Star.plots(S_exact);
    title('PosLin Exact Reach');

    subplot(2, 2, 3);
    Star.plots(S_approx);
    title('PosLin Approx Reach');

    subplot(2, 2, 4);
    Zono.plot(Z);
    title('PosLin Zono Reach');

    save_test_figure(fig, 'test_funcs_poslin', 'reach', 1, 'subdir', 'nn/poslin');

    % Save regression data
    data = struct();
    data.x_input = x;
    data.y_output = y;
    data.I_star_dim = I_star.dim;
    data.I_zono_c = I_zono.c;
    save_test_data(data, 'test_funcs_poslin', 'results', 'subdir', 'nn');
end
