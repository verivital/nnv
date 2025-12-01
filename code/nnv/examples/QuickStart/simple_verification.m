%% simple_verification.m - Your First Neural Network Verification
% This example demonstrates the core NNV workflow:
% 1. Create a neural network
% 2. Define an input set (uncertain inputs)
% 3. Compute the reachable output set
% 4. Check a safety property
%
% Usage:
%   simple_verification

fprintf('\n=== Simple Neural Network Verification ===\n\n');

%% Step 1: Create a Simple Neural Network
% A 2-input, 2-output network with one hidden layer
fprintf('Step 1: Creating neural network...\n');

% Hidden layer: 2 inputs -> 4 neurons
W1 = [1 -1; -1 1; 0.5 0.5; -0.5 -0.5];
b1 = [0; 0; 0; 0];

% Output layer: 4 neurons -> 2 outputs
W2 = [1 1 0 0; 0 0 1 1];
b2 = [0; 0];

% Build the network
layers = {
    FullyConnectedLayer(W1, b1)
    ReluLayer()
    FullyConnectedLayer(W2, b2)
};
net = NN(layers);

fprintf('  Network structure: 2 -> 4 (ReLU) -> 2\n\n');

%% Step 2: Define the Input Set
% We'll verify the network over a box of inputs
fprintf('Step 2: Defining input set...\n');

% Input bounds: x1 in [0, 1], x2 in [0, 1]
lb = [0; 0];    % lower bounds
ub = [1; 1];    % upper bounds

% Create a Star set from the box
inputBox = Box(lb, ub);
inputStar = inputBox.toStar();

fprintf('  Input: x1 in [0, 1], x2 in [0, 1]\n');
fprintf('  Input Star has %d predicate variables\n\n', inputStar.nVar);

%% Step 3: Compute Reachable Output Set
% Use approximate star reachability (fast and sound)
fprintf('Step 3: Computing reachable outputs...\n');

% Set up reachability options
reachOptions = struct;
reachOptions.reachMethod = 'approx-star';

tic;
outputStars = net.reach(inputStar, reachOptions);
reachTime = toc;

fprintf('  Reachability completed in %.3f seconds\n', reachTime);
fprintf('  Number of output sets: %d\n\n', length(outputStars));

%% Step 4: Analyze the Output
fprintf('Step 4: Analyzing output bounds...\n');

% Get bounds on each output dimension
outputBounds = outputStars.getBox();
y1_lb = outputBounds.lb(1);
y1_ub = outputBounds.ub(1);
y2_lb = outputBounds.lb(2);
y2_ub = outputBounds.ub(2);

fprintf('  Output y1 in [%.4f, %.4f]\n', y1_lb, y1_ub);
fprintf('  Output y2 in [%.4f, %.4f]\n\n', y2_lb, y2_ub);

%% Step 5: Check a Safety Property
% Safety: y1 should always be >= -1
fprintf('Step 5: Checking safety property (y1 >= -1)...\n');

safety_threshold = -1;
if y1_lb >= safety_threshold
    fprintf('  [VERIFIED] Property holds: y1 >= %.1f for all inputs\n', safety_threshold);
else
    fprintf('  [VIOLATED] Property may not hold\n');
end

%% Visualization (if 2D)
fprintf('\nStep 6: Visualization...\n');
try
    figure;

    % Plot input set
    subplot(1, 2, 1);
    inputStar.plot();
    title('Input Set');
    xlabel('x_1');
    ylabel('x_2');
    grid on;

    % Plot output set
    subplot(1, 2, 2);
    outputStars.plot();
    title('Output Reachable Set');
    xlabel('y_1');
    ylabel('y_2');
    grid on;

    fprintf('  Figure generated showing input and output sets\n');
catch
    fprintf('  (Visualization skipped - display not available)\n');
end

%% Summary
fprintf('\n=== Verification Complete ===\n');
fprintf('This example demonstrated:\n');
fprintf('  - Creating a neural network with NN()\n');
fprintf('  - Defining input sets with Box and Star\n');
fprintf('  - Computing reachable sets with net.reach()\n');
fprintf('  - Checking safety properties on outputs\n');
fprintf('\nNext: Try the MNIST example in Tutorial/NN/MNIST/\n\n');
