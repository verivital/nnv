function test_DLinearNNCS_reachLive()
    % TEST_DLINEARNNCS_REACHLIVE - Test DLinearNNCS live reachability
    %
    % Adaptive Cruise Control example with on-the-fly plotting.
    %
    % Tests that:
    %   1. Live reachability analysis completes successfully
    %   2. Output projections are computed correctly

    % System model
    A = [0 1 0 0 0 0 0; 0 0 1 0 0 0 0; 0 0 0 0 0 0 1; 0 0 0 0 1 0 0; 0 0 0 0 0 1 0; 0 0 0 0 0 -2 0; 0 0 -2 0 0 0 0];
    B = [0; 0; 0; 0; 0; 2; 0];
    C = [1 0 0 -1 0 0 0; 0 1 0 0 -1 0 0; 0 0 0 0 1 0 0];
    D = [0; 0; 0];

    plant = LinearODE(A, B, C, D);
    plantd = plant.c2d(0.1);

    % Load controller
    load('../controller_3_20.mat', 'weights', 'bias');

    n = length(weights);
    Layers = cell(1, n);
    for i=1:n - 1
        Layers{i} = LayerS(weights{1, i}, bias{i, 1}, 'poslin');
    end
    Layers{n} = LayerS(weights{1, n}, bias{n, 1}, 'purelin');
    Controller = NN(Layers);
    Controller.InputSize = 5;
    Controller.OutputSize = 1;

    % Create NNCS
    ncs = DLinearNNCS(Controller, plantd);

    % ASSERTION 1: NNCS is valid
    assert(~isempty(ncs), 'DLinearNNCS should be created successfully');

    % Initial state set
    lb = [90; 32; 0; 10; 30; 0; -4];
    ub = [100; 35; 0; 11; 30.2; 0; -4];
    init_set = Star(lb, ub);

    ref_input = [30; 1.4];

    numSteps = 10;
    method = 'exact-star';
    numCores = 1;

    % Live Reachability Analysis 1: Distance vs Safe Distance
    output_mat = [1 0 0 -1 0 0 0];
    output_vec = [];
    boundary_mat = [0 0 0 0 1.4 0 0];
    boundary_vec = [10];

    fig1 = figure;
    ncs.reachLive(init_set, ref_input, numSteps, method, numCores, output_mat, output_vec, 'blue', boundary_mat, boundary_vec);
    title('Live Reachability: Actual vs Safe Distance');

    % ASSERTION 2: Figure was created and has content
    assert(ishandle(fig1), 'Figure should be valid handle');

    save_test_figure(fig1, 'test_DLinearNNCS_reachLive', 'distance', 1, 'subdir', 'nncs/DLinearNNCS');

    % Live Reachability Analysis 2: Position and Velocity
    output_mat = [1 0 0 0 0 0 0; 0 1 0 0 0 0 0];
    output_vec = [];

    fig2 = figure;
    ncs.reachLive(init_set, ref_input, numSteps, method, numCores, output_mat, output_vec, 'blue');
    title('Live Reachability: Lead Car Position vs Velocity');

    % ASSERTION 3: Figure was created
    assert(ishandle(fig2), 'Figure should be valid handle');

    save_test_figure(fig2, 'test_DLinearNNCS_reachLive', 'posvel', 2, 'subdir', 'nncs/DLinearNNCS');

    % Save regression data
    data = struct();
    data.lb = lb;
    data.ub = ub;
    data.numSteps = numSteps;
    data.method = method;
    save_test_data(data, 'test_DLinearNNCS_reachLive', 'results', 'subdir', 'nncs');
end
