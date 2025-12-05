function test_LinearODE_sim()
    % TEST_LINEARODE_SIM - Test LinearODE point simulation methods
    %
    % Tests that:
    %   1. simDirect, simKrylov, and simOde45 produce valid results
    %   2. Different methods produce consistent results

    A = [0 1;-5 -2];
    x0 = [1; 2];
    h = 0.01;
    N = 5;

    X1 = LinearODE.simDirect(A, x0, h, N);
    X2 = LinearODE.simKrylov(A, x0, h, N, 2);
    X3 = LinearODE.simOde45(A, x0, h, N);

    % ASSERTION 1: All methods produce non-empty results
    assert(~isempty(X1), 'simDirect should produce non-empty result');
    assert(~isempty(X2), 'simKrylov should produce non-empty result');
    assert(~isempty(X3), 'simOde45 should produce non-empty result');

    % ASSERTION 2: All methods produce correct number of columns (time steps)
    assert(size(X1, 2) == N + 1, 'simDirect should produce N+1 time steps');
    assert(size(X2, 2) == N + 1, 'simKrylov should produce N+1 time steps');
    assert(size(X3, 2) == N + 1, 'simOde45 should produce N+1 time steps');

    % ASSERTION 3: All methods produce correct state dimension
    assert(size(X1, 1) == 2, 'simDirect should produce 2D state');
    assert(size(X2, 1) == 2, 'simKrylov should produce 2D state');
    assert(size(X3, 1) == 2, 'simOde45 should produce 2D state');

    % ASSERTION 4: Initial conditions match
    tol = 1e-10;
    assert(all(abs(X1(:, 1) - x0) < tol), 'simDirect initial state should match x0');
    assert(all(abs(X2(:, 1) - x0) < tol), 'simKrylov initial state should match x0');
    assert(all(abs(X3(:, 1) - x0) < tol), 'simOde45 initial state should match x0');

    % ASSERTION 5: Methods should produce similar final states
    tol_methods = 0.01;
    assert(all(abs(X1(:, end) - X3(:, end)) < tol_methods), ...
        'simDirect and simOde45 should produce similar final states');

    % Create visualization
    fig = figure;
    t = (0:N) * h;

    subplot(2, 1, 1);
    plot(t, X1(1, :), 'b-', 'LineWidth', 1.5);
    hold on;
    plot(t, X2(1, :), 'r--', 'LineWidth', 1.5);
    plot(t, X3(1, :), 'g:', 'LineWidth', 2);
    xlabel('Time (s)');
    ylabel('x_1');
    legend('Direct', 'Krylov', 'ODE45');
    title('LinearODE Simulation - State x_1');

    subplot(2, 1, 2);
    plot(t, X1(2, :), 'b-', 'LineWidth', 1.5);
    hold on;
    plot(t, X2(2, :), 'r--', 'LineWidth', 1.5);
    plot(t, X3(2, :), 'g:', 'LineWidth', 2);
    xlabel('Time (s)');
    ylabel('x_2');
    legend('Direct', 'Krylov', 'ODE45');
    title('LinearODE Simulation - State x_2');

    save_test_figure(fig, 'test_LinearODE_sim', 'comparison', 1, 'subdir', 'nncs/LinearODE');

    % Save regression data
    data = struct();
    data.A = A;
    data.x0 = x0;
    data.h = h;
    data.N = N;
    data.X1_final = X1(:, end);
    data.X2_final = X2(:, end);
    data.X3_final = X3(:, end);
    save_test_data(data, 'test_LinearODE_sim', 'results', 'subdir', 'nncs');
end
