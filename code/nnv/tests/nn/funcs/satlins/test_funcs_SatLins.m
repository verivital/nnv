function test_funcs_SatLins()
    % TEST_FUNCS_SATLINS - Test SatLins activation function operations
    %
    % Tests that:
    %   1. SatLins.evaluate produces correct output (saturates at -1 and 1)
    %   2. Various reach methods produce valid results
    %   3. Zonotope approximation works
    %   4. Step reach methods work

    %% Create input Star
    I0 = ExamplePoly.randVrep;
    I0.outerApprox;
    V = [0 0; 1 0; 0 1];
    I = Star(V', I0.A, I0.b, I0.Internal.lb, I0.Internal.ub);

    % ASSERTION 1: Input star is valid
    assert(~isempty(I), 'Input star should be created successfully');

    %% Test 1: SatLins Evaluate
    x = [-1.5; 0.5; 2];
    y = SatLins.evaluate(x);

    % ASSERTION 2: Evaluate produces correct output
    assert(isequal(y, [-1; 0.5; 1]), 'SatLins.evaluate should saturate at -1 and 1');

    %% Test 2: SatLins Reach exact star
    S_exact = SatLins.reach(I, 'exact-star');

    % ASSERTION 3: Exact reach produces valid result
    assert(~isempty(S_exact), 'exact-star reach should produce non-empty result');

    %% Test 3: SatLins Reach approx star
    S_approx = SatLins.reach(I, 'approx-star');

    % ASSERTION 4: Approx reach produces valid result
    assert(~isempty(S_approx), 'approx-star reach should produce non-empty result');

    %% Test 4: SatLins Reach abs dom
    S_absdom = SatLins.reach(I, 'abs-dom');

    % ASSERTION 5: Abstract domain reach produces valid result
    assert(~isempty(S_absdom), 'abs-dom reach should produce non-empty result');

    %% Test 5: SatLins Abstract Domain
    S_ad = SatLins.reach_abstract_domain(I);

    % ASSERTION 6: reach_abstract_domain produces valid result
    assert(~isempty(S_ad), 'reach_abstract_domain should produce non-empty result');

    %% Test 6: SatLins Reach (default)
    S_default = SatLins.reach(I);

    % ASSERTION 7: Default reach produces valid result
    assert(~isempty(S_default), 'default reach should produce non-empty result');

    %% Test 7: SatLins Reach Star Approx
    S_star_approx = SatLins.reach_star_approx(I);

    % ASSERTION 8: reach_star_approx produces valid result
    assert(~isempty(S_star_approx), 'reach_star_approx should produce non-empty result');

    %% Test 8: SatLins Reach Zono Approx
    lb = [-0.5; -0.5];
    ub = [0.5; 0.5];

    B = Box(lb, ub);
    I1 = B.toZono;

    A = [2 1; 1.5 -2];
    I_zono = I1.affineMap(A, []);

    Z = SatLins.reach_zono_approx(I_zono);

    % ASSERTION 9: Zono approx produces valid result
    assert(~isempty(Z), 'reach_zono_approx should produce non-empty result');

    I2 = I_zono.toStar;
    X = I2.sample(100);
    Y = SatLins.evaluate(X);
    S = SatLins.reach_star_approx(I2);

    fig1 = figure;
    subplot(1, 2, 1);
    Zono.plot(I_zono);
    title('Input Set');
    subplot(1, 2, 2);
    Zono.plot(Z);
    hold on;
    Star.plot(S);
    hold on;
    plot(Y(1, :), Y(2, :), '*b');
    title('Output: Star (inside) vs Zonotope (outside)');

    save_test_figure(fig1, 'test_funcs_SatLins', 'zono_vs_star', 1, 'subdir', 'nn/SatLins');

    %% Test 9: SatLins Step Reach
    S_step = SatLins.stepReach(I, 1);

    % ASSERTION 10: stepReach produces valid result
    assert(~isempty(S_step), 'stepReach should produce non-empty result');

    %% Test 10: SatLins Step Reach Abstract Domain
    S_step_ad = SatLins.stepReachAbstractDomain(I, 1);

    % ASSERTION 11: stepReachAbstractDomain produces valid result
    assert(~isempty(S_step_ad), 'stepReachAbstractDomain should produce non-empty result');

    %% Test 11: SatLins Step Reach Star Approx
    S_step_star = SatLins.stepReachStarApprox(I, 1);

    % ASSERTION 12: stepReachStarApprox produces valid result
    assert(~isempty(S_step_star), 'stepReachStarApprox should produce non-empty result');

    %% Test 12: SatLins Step Reach Zono Approx
    lb2 = [-0.5; -0.5];
    ub2 = [0.5; 0.5];

    B2 = Box(lb2, ub2);
    I1_2 = B2.toZono;

    A2 = [2 1; 1.5 -2];
    I_zono2 = I1_2.affineMap(A2, []);

    fig2 = figure;
    Zono.plot(I_zono2);
    title('Input Zonotope');

    save_test_figure(fig2, 'test_funcs_SatLins', 'input_zono', 2, 'subdir', 'nn/SatLins');

    Z2 = SatLins.stepReachZonoApprox(I_zono2, 1);
    I2_2 = I_zono2.toStar;
    S2 = SatLins.stepReachStarApprox(I2_2, 1);

    % ASSERTION 13: stepReachZonoApprox produces valid result
    assert(~isempty(Z2), 'stepReachZonoApprox should produce non-empty result');

    fig3 = figure;
    Zono.plot(Z2);
    hold on;
    Star.plot(S2);
    title('Step Reach: Zonotope vs Star');

    save_test_figure(fig3, 'test_funcs_SatLins', 'step_reach', 3, 'subdir', 'nn/SatLins');

    % Save regression data
    data = struct();
    data.input_dim = I.dim;
    data.eval_input = x;
    data.eval_output = y;
    if iscell(S_exact)
        data.num_exact_sets = length(S_exact);
    else
        data.num_exact_sets = 1;
    end
    save_test_data(data, 'test_funcs_SatLins', 'results', 'subdir', 'nn');
end
