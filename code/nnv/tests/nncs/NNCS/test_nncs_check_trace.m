function test_nncs_check_trace()
    % TEST_NNCS_CHECK_TRACE - Test NNCS.check_trace() method
    %
    % Tests that:
    %   1. NNCS can simulate a trace
    %   2. Trace checking against unsafe region works correctly

    % Load controller
    load('../controller_3_20.mat', 'weights', 'bias');
    n = length(weights);
    Layers = [];
    for i=1:n - 1
        L = LayerS(weights{1, i}, bias{i, 1}, 'poslin');
        Layers{i} = L;
    end
    Layers{n} = LayerS(weights{1, n}, bias{n, 1}, 'purelin');

    Ts = 0.1;
    reachTimeStep = 0.01;
    controlPeriod = 0.1;
    output_mat = [0 0 0 0 1 0;1 0 0 -1 0 0; 0 1 0 0 -1 0]; % feedback: relative distance, relative velocity and ego-car velocity
    Controller = NN(Layers); % feedforward neural network controller
    Controller.InputSize = 5;
    Controller.OutputSize = 1;
    Plant = NonLinearODE(6, 1, @car_dynamics, reachTimeStep, controlPeriod, output_mat);

    ncs = NNCS(Controller, Plant); % the neural network control system

    x0 = [98; 32; 0; 10; 30; 0];
    ref_input = [30; 1.4];

    N = 50;
    [simTrace, controlTrace] = ncs.evaluate(Ts, N, x0, ref_input);

    % ASSERTION 1: Simulation produces valid trace
    assert(~isempty(simTrace), 'Simulation trace should not be empty');
    assert(size(simTrace, 2) == N + 1, 'Simulation trace should have N+1 time steps');
    assert(size(simTrace, 1) == 6, 'Simulation trace should have 6 state dimensions');

    % unsafe region: dis = x1 - x2 <= 78
    unsafe_mat = [1 0 0 -1 0 0];
    unsafe_vec = [78];

    violate = NNCS.check_trace(simTrace, unsafe_mat, unsafe_vec);

    % ASSERTION 2: Check trace result is valid
    assert(violate == 0 || violate == 1, 'Violation should be 0 or 1');

    % Compute distance over time
    dis = [1 0 0 -1 0 0] * simTrace;
    times = 0:Ts:N*Ts;

    % ASSERTION 3: Distance should stay above unsafe threshold for safe trace
    min_distance = min(dis);
    if ~violate
        assert(min_distance > 78, 'Non-violating trace should have distance > 78');
    end

    % Create visualization
    fig = figure;
    plot(times, dis, '-b', 'LineWidth', 1.5);
    hold on;
    yline(78, 'r--', 'LineWidth', 1.5);
    xlabel('Time (s)');
    ylabel('Distance');
    legend('Actual Distance', 'Unsafe Threshold');
    title(sprintf('NNCS Trace Check (Violation: %d)', violate));

    save_test_figure(fig, 'test_nncs_check_trace', 'trace', 1, 'subdir', 'nncs');

    % Save regression data
    data = struct();
    data.x0 = x0;
    data.simTrace = simTrace;
    data.controlTrace = controlTrace;
    data.violate = violate;
    data.min_distance = min_distance;
    data.times = times;
    data.dis = dis;
    save_test_data(data, 'test_nncs_check_trace', 'results', 'subdir', 'nncs');
end
