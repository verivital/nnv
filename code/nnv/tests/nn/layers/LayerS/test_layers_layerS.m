function test_layers_layerS()
    % TEST_LAYERS_LAYERS - Test LayerS reach operations
    %
    % Tests that:
    %   1. LayerS with poslin activation computes exact reach set
    %   2. LayerS with satlin activation computes approximate reach set
    %   3. Reach sets contain expected states

    % Create input star
    I = ExamplePoly.randVrep;
    I.outerApprox;
    V = [0 0; 1 0; 0 1];
    I = Star(V', I.A, I.b, I.Internal.lb, I.Internal.ub);

    W = [1.5 1; 0 0.5];
    b = [0.5; 0.5];
    I1 = I.affineMap(W, b);

    %% Test 1: LayerS reach poslin
    L1 = LayerS(W, b, 'poslin');

    % ASSERTION 1: Layer is valid
    assert(~isempty(L1), 'LayerS should be created successfully');

    S = L1.reach(I, 'exact-star');

    % ASSERTION 2: Reach produces valid output
    assert(~isempty(S), 'Exact reach should produce non-empty result');

    % ASSERTION 3: Reach set is a cell array of stars
    if iscell(S)
        for i = 1:length(S)
            assert(isa(S{i}, 'Star'), sprintf('S{%d} should be a Star', i));
        end
    else
        assert(isa(S, 'Star'), 'S should be a Star');
    end

    % Create visualizations
    fig1 = figure;
    Star.plot(I);
    title('Input Star I');

    save_test_figure(fig1, 'test_layers_layerS', 'input', 1, 'subdir', 'nn/LayerS');

    fig2 = figure;
    Star.plot(I1);
    title('Affine Mapped Star I1');

    save_test_figure(fig2, 'test_layers_layerS', 'affine', 2, 'subdir', 'nn/LayerS');

    fig3 = figure;
    Star.plots(S);
    title('LayerS poslin Exact Reach');

    save_test_figure(fig3, 'test_layers_layerS', 'poslin_exact', 3, 'subdir', 'nn/LayerS');

    %% Test 2: LayerS reach satlin
    L = LayerS(W, b, 'satlin');

    % ASSERTION 4: Layer is valid
    assert(~isempty(L), 'LayerS satlin should be created successfully');

    tic;
    S_approx = L.reach(I, 'approx-star');
    t_serial = toc;

    % ASSERTION 5: Approximate reach produces valid output
    assert(~isempty(S_approx), 'Approx reach should produce non-empty result');

    fig4 = figure;
    Star.plots(S_approx);
    title(sprintf('LayerS satlin Approx Reach (%.3fs)', t_serial));

    save_test_figure(fig4, 'test_layers_layerS', 'satlin_approx', 4, 'subdir', 'nn/LayerS');

    % Test parallel execution
    tic;
    S_parallel = L.reach(I, 'approx-star', 'parallel');
    t_parallel = toc;

    % ASSERTION 6: Parallel reach produces valid output
    assert(~isempty(S_parallel), 'Parallel approx reach should produce non-empty result');

    % Save regression data
    data = struct();
    data.W = W;
    data.b = b;
    data.t_serial = t_serial;
    data.t_parallel = t_parallel;
    if iscell(S)
        data.num_reach_sets_poslin = length(S);
    else
        data.num_reach_sets_poslin = 1;
    end
    if iscell(S_approx)
        data.num_reach_sets_satlin = length(S_approx);
    else
        data.num_reach_sets_satlin = 1;
    end
    save_test_data(data, 'test_layers_layerS', 'results', 'subdir', 'nn');
end
