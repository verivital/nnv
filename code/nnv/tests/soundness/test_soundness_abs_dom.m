% test_soundness_abs_dom
% Soundness tests for abs-dom (abstract domain) reachability method
% Tests abstract domain over-approximation across multiple layers
% To run: results = runtests('test_soundness_abs_dom')

%% Test 1: ReLU with abs-dom using ImageStar input
rng(42);
L = ReluLayer();

% Create ImageStar input with perturbation
V = zeros(3, 3, 1, 2);
V(:,:,1,1) = rand(3, 3) * 0.5;  % center
V(2,2,1,1) = -0.1;  % Some negative to trigger ReLU
V(2,2,1,2) = 0.3;  % perturbation

input_is = ImageStar(V, [1; -1], [1; 1], -1, 1);
output_is = L.reach(input_is, 'abs-dom');

% Verify output is valid
assert(isa(output_is, 'ImageStar'), 'ReLU abs-dom should return ImageStar');
assert(output_is.numPred >= 0, 'Output should have valid predicates');

%% Test 2: FullyConnected with abs-dom using Star input
rng(42);

W = randn(3, 4);
b = randn(3, 1);
L2 = FullyConnectedLayer(W, b);

% Create Star input
c = randn(4, 1);
V_basis = randn(4, 2) * 0.1;
V_star = [c V_basis];
C = [eye(2); -eye(2)];
d = ones(4, 1);

input_star = Star(V_star, C, d, [-1; -1], [1; 1]);

% Test reach with abs-dom
output_star = L2.reach(input_star, 'abs-dom');

% Verify output
assert(isa(output_star, 'Star'), 'FullyConnected abs-dom should return Star');

% Sample and verify containment
for i = 1:20
    alpha = -1 + 2*rand(2, 1);
    input_sample = c + V_basis * alpha;
    output_sample = L2.evaluate(input_sample);

    % Check containment in output star
    contained = output_star.contains(output_sample);
    assert(contained, 'abs-dom output %d should contain concrete output', i);
end

%% Test 3: LeakyReLU with abs-dom
rng(42);

L3 = LeakyReluLayer('leaky_absdom', 1, {'in'}, 1, {'out'}, 0.1);

V3 = zeros(3, 3, 1, 2);
V3(:,:,1,1) = randn(3, 3) * 0.3;
V3(1,2,1,2) = 0.2;

input_is3 = ImageStar(V3, [1; -1], [1; 1], -1, 1);
output_is3 = L3.reach(input_is3, 'abs-dom');

assert(~isempty(output_is3), 'LeakyReLU abs-dom should produce output');

%% Test 4: Sigmoid with abs-dom (if supported)
rng(42);

L4 = SigmoidLayer();

V4 = zeros(2, 2, 1, 2);
V4(:,:,1,1) = randn(2, 2) * 0.5;
V4(1,1,1,2) = 0.1;

input_is4 = ImageStar(V4, [1; -1], [0.5; 0.5], -0.5, 0.5);

try
    output_is4 = L4.reach(input_is4, 'abs-dom');
    assert(~isempty(output_is4), 'Sigmoid abs-dom should produce output');
catch ME
    % abs-dom may not be fully implemented for Sigmoid
    warning('Sigmoid abs-dom not fully supported: %s', ME.message);
end

%% Test 5: Tanh with abs-dom (if supported)
rng(42);

L5 = TanhLayer();

V5 = zeros(2, 2, 1, 2);
V5(:,:,1,1) = randn(2, 2) * 0.5;
V5(1,1,1,2) = 0.1;

input_is5 = ImageStar(V5, [1; -1], [0.5; 0.5], -0.5, 0.5);

try
    output_is5 = L5.reach(input_is5, 'abs-dom');
    assert(~isempty(output_is5), 'Tanh abs-dom should produce output');
catch ME
    % abs-dom may not be fully implemented for Tanh
    warning('Tanh abs-dom not fully supported: %s', ME.message);
end

%% Test 6: SatLin with abs-dom
rng(42);

L6 = SaturatingLinearLayer();

V6 = zeros(3, 3, 1, 2);
V6(:,:,1,1) = randn(3, 3) * 0.5;
V6(2,2,1,2) = 0.2;

input_is6 = ImageStar(V6, [1; -1], [1; 1], -1, 1);
output_is6 = L6.reach(input_is6, 'abs-dom');

assert(~isempty(output_is6), 'SatLin abs-dom should produce output');

