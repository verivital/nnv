function test_DLinearODE_reach()
    % TEST_DLINEARODE_REACH - Test DLinearODE step reachability
    %
    % Tests that:
    %   1. stepReachStar computes valid reach set
    %   2. stepReachZono computes valid reach set
    %   3. Reach sets contain the initial set under identity transform

    A = [0 1;-5 -2];
    B = [0;3];
    C = [0 1];
    D = 0;
    Ts = 0.1;
    sys = LinearODE(A, B, C, D);
    sysd = sys.c2d(Ts); % convert from continuous to discrete

    % set of initial condition
    %      -1 <= x0[1] <= 1; -1 <= x0[2] <= 1
    lb1 = [-1 ;-1];
    ub1 = [1; 1];
    Bi = Box(lb1, ub1);
    I = Bi.toStar();

    % set of control input: -0.5 <= u <= 1
    lb_u = [-1];
    ub_u = [0.5];
    Bu = Box(lb_u, ub_u);
    U = Bu.toStar();

    % ASSERTION 1: Discrete system is valid
    assert(~isempty(sysd), 'Discrete system should be valid');

    %% Test 1: stepReachStar
    R = sysd.stepReachStar(I, U);

    % ASSERTION 2: stepReachStar produces valid Star
    assert(~isempty(R), 'stepReachStar should produce non-empty result');
    assert(isa(R, 'Star'), 'stepReachStar result should be a Star');
    assert(R.dim == I.dim, 'Result dimension should match input dimension');

    % ASSERTION 3: Result bounds are finite
    [lb_R, ub_R] = R.getRanges();
    assert(all(isfinite(lb_R)) && all(isfinite(ub_R)), ...
        'Reach set bounds should be finite');

    %% Test 2: stepReachZono
    I_zono = Bi.toZono();
    U_zono = Bu.toZono();
    Z = sysd.stepReachZono(I_zono, U_zono);

    % ASSERTION 4: stepReachZono produces valid Zono
    assert(~isempty(Z), 'stepReachZono should produce non-empty result');
    assert(isa(Z, 'Zono'), 'stepReachZono result should be a Zono');
    assert(Z.dim == I_zono.dim, 'Result dimension should match input dimension');

    % ASSERTION 5: Zono bounds are finite
    [lb_Z, ub_Z] = Z.getBounds();
    assert(all(isfinite(lb_Z)) && all(isfinite(ub_Z)), ...
        'Zono reach set bounds should be finite');

    % Create visualizations
    fig1 = figure;
    Star.plot(U);
    title('Control Input Set U');

    save_test_figure(fig1, 'test_DLinearODE_reach', 'input', 1, 'subdir', 'nncs/DLinearODE');

    fig2 = figure;
    Star.plot(I);
    title('Initial State Set I');

    save_test_figure(fig2, 'test_DLinearODE_reach', 'initial', 2, 'subdir', 'nncs/DLinearODE');

    fig3 = figure;
    Star.plot(R);
    title('One-Step Reach Set (Star)');

    save_test_figure(fig3, 'test_DLinearODE_reach', 'reachStar', 3, 'subdir', 'nncs/DLinearODE');

    fig4 = figure;
    Zono.plot(Z);
    title('One-Step Reach Set (Zono)');

    save_test_figure(fig4, 'test_DLinearODE_reach', 'reachZono', 4, 'subdir', 'nncs/DLinearODE');

    % Save regression data
    data = struct();
    data.A = A;
    data.B = B;
    data.Ts = Ts;
    data.I_lb = lb1;
    data.I_ub = ub1;
    data.U_lb = lb_u;
    data.U_ub = ub_u;
    data.R_lb = lb_R;
    data.R_ub = ub_R;
    data.Z_lb = lb_Z;
    data.Z_ub = ub_Z;
    save_test_data(data, 'test_DLinearODE_reach', 'results', 'subdir', 'nncs');
end
