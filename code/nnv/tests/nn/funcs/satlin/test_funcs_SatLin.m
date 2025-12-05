function test_funcs_SatLin()
    % TEST_FUNCS_SATLIN - Test SatLin activation function operations
    %
    % Tests that:
    %   1. SatLin.evaluate produces correct output (saturates at 0 and 1)
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

    %% Test 1: SatLin Evaluate
    x = [-1.5; 0.5; 2];
    y = SatLin.evaluate(x);

    % ASSERTION 2: Evaluate produces correct output
    assert(isequal(y, [0; 0.5; 1]), 'SatLin.evaluate should saturate at 0 and 1');

    %% Test 2: SatLin Reach exact star
    S_exact = SatLin.reach(I, 'exact-star');

    % ASSERTION 3: Exact reach produces valid result
    assert(~isempty(S_exact), 'exact-star reach should produce non-empty result');

    %% Test 3: SatLin Reach approx star
    S_approx = SatLin.reach(I, 'approx-star');

    % ASSERTION 4: Approx reach produces valid result
    assert(~isempty(S_approx), 'approx-star reach should produce non-empty result');

    %% Test 4: SatLin Reach abs dom
    S_absdom = SatLin.reach(I, 'abs-dom');

    % ASSERTION 5: Abstract domain reach produces valid result
    assert(~isempty(S_absdom), 'abs-dom reach should produce non-empty result');

    %% Test 5: SatLin Abstract Domain
    S_ad = SatLin.reach_abstract_domain(I);

    % ASSERTION 6: reach_abstract_domain produces valid result
    assert(~isempty(S_ad), 'reach_abstract_domain should produce non-empty result');

    %% Test 6: SatLin Reach (default)
    S_default = SatLin.reach(I);

    % ASSERTION 7: Default reach produces valid result
    assert(~isempty(S_default), 'default reach should produce non-empty result');

    %% Test 7: SatLin Reach Star Approx
    S_star_approx = SatLin.reach_star_approx(I);

    % ASSERTION 8: reach_star_approx produces valid result
    assert(~isempty(S_star_approx), 'reach_star_approx should produce non-empty result');

    %% Test 8: SatLin Reach Zono Approx
    lb = [-0.5; -0.5];
    ub = [0.5; 0.5];

    B = Box(lb, ub);
    I1 = B.toZono;

    A = [2 1; 1.5 -2];
    I_zono = I1.affineMap(A, []);

    Z = SatLin.reach_zono_approx(I_zono);

    % ASSERTION 9: Zono approx produces valid result
    assert(~isempty(Z), 'reach_zono_approx should produce non-empty result');

    I2 = I_zono.toStar;
    X = I2.sample(100);
    Y = SatLin.evaluate(X);
    S = SatLin.reach_star_approx(I2);

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

    save_test_figure(fig1, 'test_funcs_SatLin', 'zono_vs_star', 1, 'subdir', 'nn/SatLin');

    %% Test 9: SatLin Step Reach
    S_step = SatLin.stepReach(I, 1);

    % ASSERTION 10: stepReach produces valid result
    assert(~isempty(S_step), 'stepReach should produce non-empty result');

    %% Test 10: SatLin Step Reach Zono Approx
    lb2 = [-0.5; -0.5];
    ub2 = [0.5; 0.5];

    B2 = Box(lb2, ub2);
    I1_2 = B2.toZono;

    A2 = [2 1; 1.5 -2];
    I_zono2 = I1_2.affineMap(A2, []);

    fig2 = figure;
    Zono.plot(I_zono2);
    title('Input Zonotope');

    save_test_figure(fig2, 'test_funcs_SatLin', 'input_zono', 2, 'subdir', 'nn/SatLin');

    Z2 = SatLin.stepReachZonoApprox(I_zono2, 1);
    I2_2 = I_zono2.toStar;
    S2 = SatLin.stepReach(I2_2, 1);

    % ASSERTION 11: stepReachZonoApprox produces valid result
    assert(~isempty(Z2), 'stepReachZonoApprox should produce non-empty result');

    fig3 = figure;
    Zono.plot(Z2);
    hold on;
    Star.plots(S2);
    title('Step Reach: Zonotope vs Star');

    save_test_figure(fig3, 'test_funcs_SatLin', 'step_reach', 3, 'subdir', 'nn/SatLin');

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
    save_test_data(data, 'test_funcs_SatLin', 'results', 'subdir', 'nn');
end
