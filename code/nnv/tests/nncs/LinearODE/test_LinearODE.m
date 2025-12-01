function test_LinearODE()
    % TEST_LINEARODE - Test LinearODE basic simulation
    %
    % Tests that:
    %   1. LinearODE constructor works correctly
    %   2. Simulation produces valid output

    A = [0 1;-5 -2];
    B = [0;3];
    C = [0 1];
    D = 0;
    sys = LinearODE(A, B, C, D);

    % ASSERTION 1: System is valid
    assert(~isempty(sys), 'LinearODE should be created successfully');

    [u, t] = gensig('square', 4, 10, 0.1);
    x0 = [1; 2];

    % Run simulation
    [y, t_out, x] = sys.simulate(u, t, x0);

    % ASSERTION 2: Simulation produces valid output
    assert(~isempty(y), 'Simulation output y should not be empty');
    assert(~isempty(t_out), 'Simulation time t_out should not be empty');
    assert(~isempty(x), 'Simulation state x should not be empty');

    % ASSERTION 3: Output dimensions are correct
    assert(size(x, 2) == 2, 'State should have 2 dimensions');
    assert(length(t_out) == length(t), 'Output time should match input time length');

    % ASSERTION 4: State values are finite
    assert(all(isfinite(x(:))), 'All state values should be finite');

    % Create visualization
    fig = figure;
    subplot(2, 1, 1);
    plot(t_out, x(:, 1), 'b-', 'LineWidth', 1.5);
    hold on;
    plot(t_out, x(:, 2), 'r-', 'LineWidth', 1.5);
    xlabel('Time (s)');
    ylabel('State');
    legend('x_1', 'x_2');
    title('LinearODE State Trajectories');

    subplot(2, 1, 2);
    plot(t_out, y, 'g-', 'LineWidth', 1.5);
    hold on;
    plot(t, u, 'k--', 'LineWidth', 1);
    xlabel('Time (s)');
    ylabel('Output / Input');
    legend('Output y', 'Input u');
    title('LinearODE Output');

    save_test_figure(fig, 'test_LinearODE', 'simulate', 1, 'subdir', 'nncs/LinearODE');

    % Save regression data
    data = struct();
    data.A = A;
    data.B = B;
    data.C = C;
    data.D = D;
    data.x0 = x0;
    data.x_final = x(end, :)';
    data.y_final = y(end);
    save_test_data(data, 'test_LinearODE', 'results', 'subdir', 'nncs');
end
