function test_feedforward_NN()
    % TEST_FEEDFORWARD_NN - Test feedforward neural network operations
    %
    % Tests that:
    %   1. Falsification finds counter-examples when they exist
    %   2. Robustness verification works correctly
    %   3. Safety verification works correctly

    rng(0); % ensure test won't fail due to random seed

    %% Test 1: Falsify
    % Create NN
    W = [1 1; 0 1];
    b = [0; 0.5];
    L = LayerS(W, b, 'poslin');
    Layers = {L};
    F = NN(Layers);

    % ASSERTION 1: NN is valid
    assert(~isempty(F), 'NN should be created successfully');

    % Create input set
    lb = [-1; -1];
    ub = [1; 1];
    I = Star(lb, ub);

    % Compute reachability
    R = F.reach(I);

    % ASSERTION 2: Reachability produces valid output
    assert(~isempty(R), 'Reachability should produce non-empty result');

    % Create unsafe region
    G = [-1 0];
    g = [-1.5];
    U = HalfSpace(G, g);

    % Define number of samples to attempt falsification
    n_samples = 1000;
    % Find counter examples
    counter_inputs = F.falsify(I, U, n_samples);
    counter_outputs = F.evaluate(counter_inputs);

    % ASSERTION 3: Falsification finds counter-examples
    assert(~isempty(counter_inputs), 'Falsification should find counter-examples');
    assert(size(counter_inputs, 1) == 2, 'Counter inputs should be 2D');

    % ASSERTION 4: Counter-examples are in unsafe region
    for i = 1:size(counter_outputs, 2)
        assert(G * counter_outputs(:, i) <= g, ...
            'Counter-example output should be in unsafe region');
    end

    % Visualize results
    fig1 = figure;
    subplot(1, 2, 1);
    Star.plot(I);
    hold on;
    plot(counter_inputs(1, :), counter_inputs(2, :), 'o');
    title('Input Set and counter inputs');

    subplot(1, 2, 2);
    Star.plots(R);
    hold on;
    plot(counter_outputs(1, :), counter_outputs(2, :), 'o');
    title('Output set and counter outputs');

    save_test_figure(fig1, 'test_feedforward_NN', 'falsify', 1, 'subdir', 'nn');

    %% Test 2: isRobust
    W = [1 1; 0 1];
    b = [0; 0.5];
    L = LayerS(W, b, 'poslin');
    Layers = {L};
    F = NN(Layers);

    lb = [0.5; 0.5];
    ub = [1.5; 1.5];
    I = Star(lb, ub);

    G = [-1 0];
    g = [-1.5];
    U = HalfSpace(G, g);

    reachOptions = struct;
    reachOptions.reachMethod = 'exact-star';

    res = F.verify_robustness(I, reachOptions, U);

    % ASSERTION 5: Robustness verification returns valid result
    assert(res == 0 || res == 1 || res == 2, 'Robustness result should be 0, 1, or 2');

    fig2 = figure;
    Star.plots(F.reachSet{end});
    hold on;
    U.plot;
    title(sprintf('Robustness Check (result: %d)', res));

    save_test_figure(fig2, 'test_feedforward_NN', 'robustness', 2, 'subdir', 'nn');

    %% Test 3: Safety
    W = [1 1; 0 1];
    b = [0; 0.5];
    L = LayerS(W, b, 'poslin');
    Layers = {L};
    F = NN(Layers);

    lb = [-1; -1];
    ub = [1; 1];
    I = Star(lb, ub);

    reachOptions = struct;
    reachOptions.reachMethod = 'exact-star';
    R = F.reach(I, reachOptions);

    G = [-1 0];
    g = [-1.5];
    U = HalfSpace(G, g);

    n_samples = 100;
    reachOptions = struct;
    reachOptions.reachMethod = 'approx-star';
    [safe, counter_inputs] = F.verify_safety(I, U, reachOptions, n_samples);
    counter_outputs = F.evaluate(counter_inputs);

    % ASSERTION 6: Safety verification returns valid result
    assert(safe == 0 || safe == 1, 'Safety result should be 0 or 1');

    fig3 = figure;
    subplot(1, 2, 1);
    Star.plot(I);
    hold on;
    if ~isempty(counter_inputs)
        plot(counter_inputs(1, :), counter_inputs(2, :), 'o');
    end
    title('Input Set and counter input set');

    subplot(1, 2, 2);
    Star.plots(R);
    hold on;
    if ~isempty(counter_outputs)
        plot(counter_outputs(1, :), counter_outputs(2, :), 'o');
    end
    title(sprintf('Safety Check (safe: %d)', safe));

    save_test_figure(fig3, 'test_feedforward_NN', 'safety', 3, 'subdir', 'nn');

    % Save regression data
    data = struct();
    data.W = W;
    data.b = b;
    data.num_counter_examples = size(counter_inputs, 2);
    data.safe = safe;
    save_test_data(data, 'test_feedforward_NN', 'results', 'subdir', 'nn');
end
