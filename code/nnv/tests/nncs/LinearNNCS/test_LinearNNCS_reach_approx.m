function test_LinearNNCS_reach_approx()
    % TEST_LINEARNNCS_REACH_APPROX - Test LinearNNCS approximate reachability
    %
    % Adaptive Cruise Control example with continuous linear plant.
    %
    % Tests that:
    %   1. Approximate reachability analysis completes successfully
    %   2. Output reach sets are computed

    % System model
    A = [0 1 0 0 0 0 0; 0 0 1 0 0 0 0; 0 0 0 0 0 0 1; 0 0 0 0 1 0 0; 0 0 0 0 0 1 0; 0 0 0 0 0 -2 0; 0 0 -2 0 0 0 0];
    B = [0; 0; 0; 0; 0; 2; 0];
    C = [1 0 0 -1 0 0 0; 0 1 0 0 -1 0 0; 0 0 0 0 1 0 0];
    D = [0; 0; 0];

    plant = LinearODE(A, B, C, D);

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
    ncs = LinearNNCS(Controller, plant);

    % ASSERTION 1: NNCS is valid
    assert(~isempty(ncs), 'LinearNNCS should be created successfully');

    % Initial state set
    lb = [90; 32; 0; 10; 30; 0; -4];
    ub = [100; 35; 0; 11; 30.2; 0; -4];

    ref_input = [30; 1.4];

    % Reachability Analysis
    reachPRM.init_set = Star(lb, ub);
    reachPRM.ref_input = ref_input;
    reachPRM.numSteps = 10;
    reachPRM.reachMethod = 'approx-star';
    reachPRM.numCores = 1;

    [R, reachTime] = ncs.reach(reachPRM);

    % ASSERTION 2: Reachability produces valid output
    assert(~isempty(R), 'Reachability should produce non-empty result');
    assert(reachTime > 0, 'Reach time should be positive');

    % Create visualizations
    fig1 = figure;
    map_mat = [1 0 0 -1 0 0 0];
    map_vec = [];
    ncs.plotOutputReachSets('blue', map_mat, map_vec);
    hold on;
    map_mat = [0 0 0 0 1.4 0 0];
    map_vec = [10];
    ncs.plotOutputReachSets('red', map_mat, map_vec);
    title('Approx: Actual Distance (blue) vs Safe Distance (red)');
    xlabel('Time Step');
    ylabel('Distance');

    save_test_figure(fig1, 'test_LinearNNCS_reach_approx', 'distance', 1, 'subdir', 'nncs/LinearNNCS');

    fig2 = figure;
    map_mat = [1 0 0 -1 0 0 0; 0 0 0 0 1 0 0];
    map_vec = [];
    ncs.plotOutputReachSets('blue', map_mat, map_vec);
    title('Approx: Distance vs Ego Car Speed');
    xlabel('Distance');
    ylabel('Speed');

    save_test_figure(fig2, 'test_LinearNNCS_reach_approx', 'dist_speed', 2, 'subdir', 'nncs/LinearNNCS');

    % Save regression data
    data = struct();
    data.lb = lb;
    data.ub = ub;
    data.numSteps = reachPRM.numSteps;
    data.reachMethod = reachPRM.reachMethod;
    data.reachTime = reachTime;
    save_test_data(data, 'test_LinearNNCS_reach_approx', 'results', 'subdir', 'nncs');
end
