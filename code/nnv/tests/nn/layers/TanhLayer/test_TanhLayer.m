% Test TanhLayer functionality
% To run: results = runtests('test_TanhLayer')

%% Test 1: TanhLayer constructor - default
L = TanhLayer();
assert(strcmp(L.Name, 'act_func_layer'));

%% Test 2: TanhLayer constructor - with name
L = TanhLayer({'tanh1'});
assert(strcmp(L.Name, 'tanh1'));

%% Test 3: TanhLayer evaluate - simple input
L = TanhLayer();

% Test with simple matrix
input = [-2 -1 0 1 2; -3 -0.5 0.5 3 4];
output = L.evaluate(input);

% Tanh should be applied element-wise
expected = tanh(input);
assert(all(abs(output(:) - expected(:)) < 1e-6), 'TanhLayer evaluate failed');

%% Test 4: TanhLayer evaluate - 3D image
L = TanhLayer();

% Create 3D image
IM(:,:,1) = [-1 1 0 -1; 0 0 1 -1; 1 0 -1 0; 1 -1 -1 1];
IM(:,:,2) = [0 1 0 0; 1 0 0 -1; 0 1 -1 0; 0 1 0 -1];
IM(:,:,3) = [1 -1 1 1; 1 -1 0 1; 0 1 -1 0; 1 0 -1 0];

output = L.evaluate(IM);

% Tanh should be applied element-wise
expected = tanh(IM);
assert(all(abs(output(:) - expected(:)) < 1e-6), 'TanhLayer 3D evaluate failed');

%% Test 5: TanhLayer evaluateSequence
L = TanhLayer();

% Create sequence input
seq_input = rand(10, 5, 3);
output_seq = L.evaluateSequence(seq_input);

% Tanh should be applied element-wise
expected_seq = tanh(seq_input);
assert(all(abs(output_seq(:) - expected_seq(:)) < 1e-6), 'TanhLayer evaluateSequence failed');

%% Test 6: TanhLayer reach with ImageStar - approx-star
L = TanhLayer();

% Create ImageStar input
IM(:,:,1) = [1 1 0 1; 0 0 1 1; 1 0 1 0; 1 1 1 1];
IM(:,:,2) = [0 1 0 0; 1 0 0 1; 0 1 1 0; 0 0 0 1];
IM(:,:,3) = [1 1 1 1; 1 1 0 1; 0 1 1 0; 1 0 1 0];

LB(:,:,1) = [-0.1 -0.1 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
LB(:,:,2) = [-0.1 -0.1 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
LB(:,:,3) = LB(:,:,2);

UB(:,:,1) = [0.1 0.1 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
UB(:,:,2) = [0.1 0.1 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
UB(:,:,3) = UB(:,:,2);

image_star = ImageStar(IM, LB, UB);
images = L.reach(image_star, 'approx-star');

assert(~isempty(images), 'TanhLayer reach approx-star should return result');
assert(isa(images(1), 'ImageStar'), 'reach should return ImageStar');

%% Test 7: TanhLayer reach with ImageZono
L = TanhLayer();

LB(:,:,1) = [-0.1 -0.1 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
LB(:,:,2) = [-0.1 -0.1 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
LB(:,:,3) = LB(:,:,2);

UB(:,:,1) = [0.1 0.1 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
UB(:,:,2) = [0.1 0.1 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
UB(:,:,3) = UB(:,:,2);

image_zono = ImageZono(LB, UB);
image = L.reach_zono(image_zono);

assert(isa(image, 'ImageZono'), 'reach_zono should return ImageZono');

%% Test 8: TanhLayer reach with Star - approx-star
L = TanhLayer();

% Create Star input
lb = [-1; -1];
ub = [1; 1];
B = Box(lb, ub);
I_star = B.toStar;

S = L.reach_star_single_input(I_star, 'approx-star');
assert(~isempty(S), 'TanhLayer reach with Star should return result');

%% Test 9: TanhLayer toGPU
L = TanhLayer();
L_gpu = L.toGPU();
assert(isa(L_gpu, 'TanhLayer'));

%% Test 10: TanhLayer changeParamsPrecision
L = TanhLayer();
L_single = L.changeParamsPrecision('single');
assert(isa(L_single, 'TanhLayer'));
