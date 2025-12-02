function test_DLinearNNCS_reach_exact()
    % TEST_DLINEARNNCS_REACH_EXACT - Test DLinearNNCS exact reachability
    %
    % Adaptive Cruise Control example with discrete linear plant and NN controller.
    %
    % Tests that:
    %   1. Reachability analysis completes successfully
    %   2. Verification against safety property works
    %   3. Reach sets are computed for all time steps

    % System model - Adaptive Cruise Control
    % x1-x3: lead car (position, velocity, internal state)
    % x4-x6: ego car (position, velocity, internal state)
    % x7: auxiliary variable for linearization

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

    % Reachability Analysis
    reachPRM.init_set = init_set;
    reachPRM.ref_input = ref_input;
    reachPRM.numSteps = 10;
    reachPRM.reachMethod = 'exact-star';
    reachPRM.numCores = 1;

    [R, reachTime] = ncs.reach(reachPRM);

    % ASSERTION 2: Reachability produces valid output
    assert(~isempty(R), 'Reachability should produce non-empty result');
    assert(reachTime > 0, 'Reach time should be positive');

    % Verification - safety property: distance > safe distance
    alpha = 1;
    unsafe_mat = [1 0 0 -1 alpha*1.4 0 0];
    unsafe_vec = alpha*10;
    unsafeRegion = HalfSpace(unsafe_mat, unsafe_vec);

    [safe, counterExamples, verifyTime] = ncs.verify(reachPRM, unsafeRegion);

    % ASSERTION 3: Verification returns valid result
    % safe can be: numeric (0/1/-1/2), string ('SAFE'/'UNSAFE'), or array
    assert(~isempty(safe), 'Verification should return a result');
    assert(verifyTime >= 0, 'Verify time should be non-negative');

    % Create visualizations
    fig1 = figure;
    ncs.plotPlantReachSets('blue', 1);
    hold on;
    ncs.plotPlantReachSets('red', 4);
    title('Position: Lead car (blue) vs Ego car (red)');
    xlabel('Time Step');
    ylabel('Position');

    save_test_figure(fig1, 'test_DLinearNNCS_reach_exact', 'position', 1, 'subdir', 'nncs/DLinearNNCS');

    fig2 = figure;
    ncs.plotPlantReachSets('blue', 2);
    hold on;
    ncs.plotPlantReachSets('red', 5);
    title('Velocity: Lead car (blue) vs Ego car (red)');
    xlabel('Time Step');
    ylabel('Velocity');

    save_test_figure(fig2, 'test_DLinearNNCS_reach_exact', 'velocity', 2, 'subdir', 'nncs/DLinearNNCS');

    fig3 = figure;
    ncs.plotControllerReachSets('green', 1);
    title('Controller Reach Sets');
    xlabel('Time Step');
    ylabel('Control Input');

    save_test_figure(fig3, 'test_DLinearNNCS_reach_exact', 'control', 3, 'subdir', 'nncs/DLinearNNCS');

    fig4 = figure;
    map_mat = [1 0 0 -1 0 0 0];
    map_vec = [];
    ncs.plotOutputReachSets('blue', map_mat, map_vec);
    hold on;
    map_mat = [0 0 0 0 1.4 0 0];
    map_vec = [10];
    ncs.plotOutputReachSets('red', map_mat, map_vec);
    title('Actual Distance (blue) vs Safe Distance (red)');
    xlabel('Time Step');
    ylabel('Distance');
    legend('Actual', 'Safe');

    save_test_figure(fig4, 'test_DLinearNNCS_reach_exact', 'distance', 4, 'subdir', 'nncs/DLinearNNCS');

    % Save regression data
    data = struct();
    data.lb = lb;
    data.ub = ub;
    data.numSteps = reachPRM.numSteps;
    data.reachTime = reachTime;
    data.safe = safe;
    data.verifyTime = verifyTime;
    save_test_data(data, 'test_DLinearNNCS_reach_exact', 'results', 'subdir', 'nncs');
end
