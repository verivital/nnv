% test_regression_neural_ode
% Regression tests for Neural ODE verification
% Tests ODEblockLayer and related ODE classes
% To run: results = runtests('test_regression_neural_ode')

%% Test 1: ODEblockLayer class exists
rng(42);
assert(exist('ODEblockLayer', 'class') == 8, 'ODEblockLayer class should exist');

%% Test 2: LinearODE class exists
rng(42);
assert(exist('LinearODE', 'class') == 8, 'LinearODE class should exist');

%% Test 3: NonLinearODE class exists
rng(42);
assert(exist('NonLinearODE', 'class') == 8, 'NonLinearODE class should exist');

%% Test 4: Create simple LinearODE model
rng(42);
A = [-1 0; 0 -2];
B = [1; 0];
C = eye(2);
D = zeros(2,1);
controlPeriod = 1;
numSteps = 10;
ode = LinearODE(A, B, C, D, controlPeriod, numSteps);
assert(~isempty(ode), 'LinearODE should be created');

%% Test 5: Create ODEblockLayer from LinearODE
rng(42);
A = [-1 0; 0 -2];
B = [1; 0];
C = eye(2);
D = zeros(2,1);
controlPeriod = 1;
numSteps = 10;
ode = LinearODE(A, B, C, D, controlPeriod, numSteps);
layer = ODEblockLayer(ode, controlPeriod, 0.1, false);
assert(~isempty(layer), 'ODEblockLayer should be created');

%% Test 6: ODEblockLayer methods via metaclass
rng(42);
mc = ?ODEblockLayer;
methods = {mc.MethodList.Name};
assert(ismember('evaluate', methods), 'ODEblockLayer should have evaluate method');
assert(ismember('reach', methods), 'ODEblockLayer should have reach method');

%% Test 7: ODEblockLayer evaluate works
rng(42);
A = [-1 0; 0 -2];
B = [1; 0];
C = eye(2);
D = zeros(2,1);
controlPeriod = 1;
numSteps = 10;
ode = LinearODE(A, B, C, D, controlPeriod, numSteps);
layer = ODEblockLayer(ode, controlPeriod, 0.1, false);
x0 = [1; 0.5];
y = layer.evaluate(x0);
assert(length(y) == 2, 'Output should have 2 elements');
assert(all(isfinite(y)), 'Output should be finite');

%% Test 8: Load spiral model if available
rng(42);
model_path = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, ...
    'NN', filesep, 'NeuralODEs', filesep, 'models', filesep, 'odeffnn_spiral.mat'];
if exist(model_path, 'file')
    net_info = load(model_path);
    assert(isfield(net_info, 'Wb'), 'Model should have Wb field');
else
    assert(true, 'Model file not found, skipping');
end
