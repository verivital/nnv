%% test_installation.m - Verify NNV Installation
% Run this script to check that NNV is properly installed and configured.
%
% Usage:
%   test_installation
%
% Expected output:
%   All tests should show [PASS]

fprintf('\n=== NNV Installation Test ===\n\n');

all_passed = true;

%% Test 1: Check NNV Version
fprintf('Test 1: NNV Version... ');
try
    v = NNVVERSION();
    fprintf('[PASS] %s\n', v);
catch
    fprintf('[FAIL] NNVVERSION not found\n');
    all_passed = false;
end

%% Test 2: Star Class Available
fprintf('Test 2: Star class... ');
if exist('Star', 'class') == 8
    fprintf('[PASS]\n');
else
    fprintf('[FAIL] Star class not on path\n');
    all_passed = false;
end

%% Test 3: Create a Star Set
fprintf('Test 3: Create Star set... ');
try
    % Create a 2D Star set from lower and upper bounds
    lb = [-1; -1];  % lower bounds
    ub = [1; 1];    % upper bounds
    S = Star(lb, ub);
    fprintf('[PASS]\n');
catch ME
    fprintf('[FAIL] %s\n', ME.message);
    all_passed = false;
end

%% Test 4: Box to Star Conversion
fprintf('Test 4: Box to Star... ');
try
    lb = [-1; -1];  % lower bounds
    ub = [1; 1];    % upper bounds
    B = Box(lb, ub);
    S = B.toStar();
    fprintf('[PASS]\n');
catch ME
    fprintf('[FAIL] %s\n', ME.message);
    all_passed = false;
end

%% Test 5: FullyConnectedLayer
fprintf('Test 5: FullyConnectedLayer... ');
try
    W = [1 2; 3 4];
    b = [0; 0];
    layer = FullyConnectedLayer(W, b);
    out = layer.evaluate([1; 1]);
    expected = [3; 7];
    if isequal(out, expected)
        fprintf('[PASS]\n');
    else
        fprintf('[FAIL] Wrong output\n');
        all_passed = false;
    end
catch ME
    fprintf('[FAIL] %s\n', ME.message);
    all_passed = false;
end

%% Test 6: ReLU Layer
fprintf('Test 6: ReLU Layer... ');
try
    relu = ReluLayer();
    out = relu.evaluate([-2; 0; 3]);
    expected = [0; 0; 3];
    if isequal(out, expected)
        fprintf('[PASS]\n');
    else
        fprintf('[FAIL] Wrong output\n');
        all_passed = false;
    end
catch ME
    fprintf('[FAIL] %s\n', ME.message);
    all_passed = false;
end

%% Test 7: Simple Network Creation
fprintf('Test 7: Create NN... ');
try
    rng(42); % set random seed for reproducibility
    W1 = rand(4, 2);
    b1 = rand(4, 1);
    W2 = rand(2, 4);
    b2 = rand(2, 1);

    layers = {
        FullyConnectedLayer(W1, b1)
        ReluLayer()
        FullyConnectedLayer(W2, b2)
    };

    net = NN(layers);
    fprintf('[PASS]\n');
catch ME
    fprintf('[FAIL] %s\n', ME.message);
    all_passed = false;
end

%% Test 8: Network Evaluation
fprintf('Test 8: NN evaluate... ');
try
    % Use the network from Test 7
    input = [0.5; 0.5];
    output = net.evaluate(input);
    if length(output) == 2
        fprintf('[PASS]\n');
    else
        fprintf('[FAIL] Wrong output dimension\n');
        all_passed = false;
    end
catch ME
    fprintf('[FAIL] %s\n', ME.message);
    all_passed = false;
end

%% Summary
fprintf('\n=== Summary ===\n');
if all_passed
    fprintf('All tests PASSED! NNV is ready to use.\n');
    fprintf('\nNext steps:\n');
    fprintf('  - Run simple_verification.m for your first verification\n');
    fprintf('  - Explore examples/Tutorial/ for more examples\n');
else
    fprintf('Some tests FAILED.\n');
    fprintf('Run check_nnv_setup() for detailed diagnostics.\n');
end
fprintf('\n');
