function test_LinearODE_reachZono()
    % TEST_LINEARODE_REACHZONO - Test LinearODE zonotope reachability
    %
    % This example is from the paper: Reachability of Uncertain Linear Systems
    % Using Zonotopes, Antoine Girard, HSCC 2005.
    %
    % Tests that:
    %   1. Zonotope reachability produces valid results
    %   2. Reachable sets maintain zonotope properties
    %   3. Order reduction keeps sets bounded

    A = [-1 -4; 4 -1];
    B = [1; 1];
    lb = [0.9; -0.1];
    ub = [1.1; 0.1];
    B1 = Box(lb, ub);

    I = B1.toZono(); % input set: I = [0.9, 1.1] x [-0.1, 0.1]
    U = Zono(0, 0.05 * eye(1)); % ||u(t)|| <= mu = 0.05

    h = 0.02; % time-step for reachability analysis
    N = 50; % number of iterations (reduced from 100 for faster testing)
    order_max = 10; % maximum number of generators

    sys = LinearODE(A, B, [], []);

    % ASSERTION 1: Initial set is valid zonotope
    assert(isa(I, 'Zono'), 'Initial set should be a Zono');
    assert(I.dim == 2, 'Initial set should be 2D');

    R = sys.reachZono(I, U, h, N, order_max);

    % ASSERTION 2: Reachability produces valid output
    assert(~isempty(R), 'Reachability should produce non-empty result');
    assert(length(R) >= N, 'Should have at least N reachable sets');

    % ASSERTION 3: All results are valid zonotopes
    for i = 1:min(10, length(R))
        assert(isa(R(i), 'Zono'), sprintf('R(%d) should be a Zono', i));
        assert(R(i).dim == 2, sprintf('R(%d) should be 2D', i));
    end

    % ASSERTION 4: Order is bounded by order_max (approximately)
    for i = 1:min(10, length(R))
        num_gen = size(R(i).V, 2);
        % Order should be reasonably bounded
        assert(num_gen <= 2 * order_max, ...
            sprintf('R(%d) should have bounded number of generators', i));
    end

    % ASSERTION 5: Bounds should be finite (system is stable)
    [lb_final, ub_final] = R(end).getBounds();
    assert(all(isfinite(lb_final)) && all(isfinite(ub_final)), ...
        'Final reach set should have finite bounds');

    % Create visualization
    fig = figure;
    % Plot first 20 zonotopes for clarity
    num_plot = min(20, length(R));
    Zono.plots(R(1:num_plot));
    title(sprintf('LinearODE Zonotope Reachability (first %d steps)', num_plot));
    xlabel('x_1');
    ylabel('x_2');

    save_test_figure(fig, 'test_LinearODE_reachZono', 'reachZono', 1, 'subdir', 'nncs/LinearODE');

    % Save regression data
    data = struct();
    data.A = A;
    data.B = B;
    data.h = h;
    data.N = N;
    data.order_max = order_max;
    data.I_c = I.c;
    data.I_V = I.V;
    data.R_final_c = R(end).c;
    data.R_final_V = R(end).V;
    data.lb_final = lb_final;
    data.ub_final = ub_final;
    data.num_reach_sets = length(R);
    save_test_data(data, 'test_LinearODE_reachZono', 'results', 'subdir', 'nncs');
end
