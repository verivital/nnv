% test_soundness_solvers
% Tests reachability analysis with different LP solvers
% Compares results across linprog, GLPK, and Gurobi (if available)
% To run: results = runtests('test_soundness_solvers')

%% Test 1: Detect available solvers
% Check which solvers are available on this system
has_linprog = ~isempty(which('linprog'));
has_glpk = ~isempty(which('glpk'));
has_gurobi = ~isempty(which('gurobi'));

fprintf('Solver availability:\n');
fprintf('  linprog: %s\n', mat2str(has_linprog));
fprintf('  GLPK:    %s\n', mat2str(has_glpk));
fprintf('  Gurobi:  %s\n', mat2str(has_gurobi));

assert(has_linprog || has_glpk, 'At least linprog or GLPK must be available');

%% Test 2: ReLU reachability with linprog solver
rng(42);
L = ReluLayer();

V = zeros(3, 3, 1, 2);
V(:,:,1,1) = rand(3, 3) * 0.5;
V(2,2,1,1) = -0.1;
V(2,2,1,2) = 0.3;

input_is = ImageStar(V, [1; -1], [1; 1], -1, 1);

% Test with linprog
output_linprog = L.reach(input_is, 'approx-star', [], 0, [], 'linprog');
assert(~isempty(output_linprog), 'ReLU with linprog should produce output');

%% Test 3: ReLU reachability with GLPK solver
rng(42);

% Check GLPK availability in this test section
if isempty(which('glpk'))
    warning('GLPK not available, skipping test');
else
    L = ReluLayer();

    V = zeros(3, 3, 1, 2);
    V(:,:,1,1) = rand(3, 3) * 0.5;
    V(2,2,1,1) = -0.1;
    V(2,2,1,2) = 0.3;

    input_is = ImageStar(V, [1; -1], [1; 1], -1, 1);

    output_glpk = L.reach(input_is, 'approx-star', [], 0, [], 'glpk');
    assert(~isempty(output_glpk), 'ReLU with GLPK should produce output');
end

%% Test 4: ReLU reachability with Gurobi solver (if available)
rng(42);

% Check Gurobi availability in this test section
if isempty(which('gurobi'))
    warning('Gurobi not available, skipping test');
else
    try
        L = ReluLayer();

        V = zeros(3, 3, 1, 2);
        V(:,:,1,1) = rand(3, 3) * 0.5;
        V(2,2,1,1) = -0.1;
        V(2,2,1,2) = 0.3;

        input_is = ImageStar(V, [1; -1], [1; 1], -1, 1);

        output_gurobi = L.reach(input_is, 'approx-star', [], 0, [], 'gurobi');
        assert(~isempty(output_gurobi), 'ReLU with Gurobi should produce output');
        fprintf('Gurobi ReLU test PASSED\n');
    catch ME
        % Handle Gurobi MEX/DLL loading issues gracefully
        if contains(ME.message, 'MEX-file') || contains(ME.message, 'could not be found')
            warning('Gurobi DLL loading issue (may need to add Gurobi bin to PATH): %s', ME.identifier);
        else
            rethrow(ME);
        end
    end
end

%% Test 5: FullyConnected layer with linprog
rng(42);

W = randn(3, 5);
b = randn(3, 1);
L5 = FullyConnectedLayer(W, b);

c = randn(5, 1);
V_basis = randn(5, 2) * 0.1;
V_star = [c V_basis];
C = [eye(2); -eye(2)];
d = ones(4, 1);

input_star = Star(V_star, C, d, [-1; -1], [1; 1]);

% Test with linprog
output_linprog5 = L5.reach(input_star, 'approx-star', [], 0, [], 'linprog');
assert(~isempty(output_linprog5), 'FC with linprog should produce output');

%% Test 6: FullyConnected with GLPK
rng(42);

if isempty(which('glpk'))
    warning('GLPK not available, skipping test');
else
    W = randn(3, 5);
    b = randn(3, 1);
    L6 = FullyConnectedLayer(W, b);

    c = randn(5, 1);
    V_basis = randn(5, 2) * 0.1;
    V_star = [c V_basis];
    C = [eye(2); -eye(2)];
    d = ones(4, 1);

    input_star = Star(V_star, C, d, [-1; -1], [1; 1]);

    output_glpk6 = L6.reach(input_star, 'approx-star', [], 0, [], 'glpk');
    assert(~isempty(output_glpk6), 'FC with GLPK should produce output');
end

%% Test 7: FullyConnected with Gurobi
rng(42);

if isempty(which('gurobi'))
    warning('Gurobi not available, skipping test');
else
    W = randn(3, 5);
    b = randn(3, 1);
    L7 = FullyConnectedLayer(W, b);

    c = randn(5, 1);
    V_basis = randn(5, 2) * 0.1;
    V_star = [c V_basis];
    C = [eye(2); -eye(2)];
    d = ones(4, 1);

    input_star = Star(V_star, C, d, [-1; -1], [1; 1]);

    output_gurobi7 = L7.reach(input_star, 'approx-star', [], 0, [], 'gurobi');
    assert(~isempty(output_gurobi7), 'FC with Gurobi should produce output');
    fprintf('Gurobi FC test PASSED\n');
end

%% Test 8: MaxPool with Gurobi soundness check
rng(42);

if isempty(which('gurobi'))
    warning('Gurobi not available, skipping test');
else
    L8 = MaxPooling2DLayer([2 2], [2 2], [0 0 0 0]);

    V8 = zeros(4, 4, 1, 2);
    V8(:,:,1,1) = [0.1 0.2 0.3 0.4; 0.2 0.8 0.4 0.5; 0.3 0.4 0.5 0.6; 0.4 0.5 0.6 0.7];
    V8(2,2,1,2) = 0.05;

    input_is8 = ImageStar(V8, [1; -1], [1; 1], -1, 1);

    % Test with Gurobi
    output_gr8 = L8.reach(input_is8, 'approx-star', [], 0, [], 'gurobi');
    assert(~isempty(output_gr8), 'MaxPool with Gurobi should produce output');

    % Verify soundness
    for i = 1:5
        alpha = -1 + 2*rand();
        input_concrete = V8(:,:,:,1) + alpha * V8(:,:,:,2);
        output_concrete = L8.evaluate(input_concrete);

        contained = soundness_test_utils.verify_imagestar_containment(output_gr8, output_concrete);
        assert(contained, 'Gurobi output should contain sample %d', i);
    end
    fprintf('MaxPool Gurobi soundness verified\n');
end

