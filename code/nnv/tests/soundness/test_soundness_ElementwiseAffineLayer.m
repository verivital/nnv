% test_soundness_ElementwiseAffineLayer
% Tests for Elementwise Affine layer
% To run: results = runtests('test_soundness_ElementwiseAffineLayer')

%% Test 1: ElementwiseAffineLayer constructor (4 args)
rng(42);

scale = 2.0;
offset = 0.5;
doScale = true;
doOffset = true;

L = ElementwiseAffineLayer(scale, offset, doScale, doOffset);

assert(L.Scale == scale, 'Scale should match');
assert(L.Offset == offset, 'Offset should match');
assert(L.DoScale == true, 'DoScale should be true');
assert(L.DoOffset == true, 'DoOffset should be true');

%% Test 2: ElementwiseAffineLayer constructor (5 args with name)
rng(42);

L2 = ElementwiseAffineLayer('affine_test', 1.5, 0.25, true, true);

assert(strcmp(L2.Name, 'affine_test'), 'Name should be set correctly');
assert(L2.Scale == 1.5, 'Scale should be 1.5');
assert(L2.Offset == 0.25, 'Offset should be 0.25');

%% Test 3: ElementwiseAffineLayer evaluate with scale only
rng(42);

scale = 2.0;
offset = 0.0;
L3 = ElementwiseAffineLayer(scale, offset, true, false);

input = rand(3, 3);
output = L3.evaluate(input);

expected = input * scale;
assert(max(abs(output(:) - expected(:))) < 1e-10, 'Scale-only should multiply input');

%% Test 4: ElementwiseAffineLayer evaluate with offset only
rng(42);

scale = 1.0;
offset = 0.5;
L4 = ElementwiseAffineLayer(scale, offset, false, true);

input = rand(3, 3);
output = L4.evaluate(input);

expected = input + offset;
assert(max(abs(output(:) - expected(:))) < 1e-10, 'Offset-only should add to input');

%% Test 5: ElementwiseAffineLayer evaluate with both
rng(42);

scale = 2.0;
offset = 0.5;
L5 = ElementwiseAffineLayer(scale, offset, true, true);

input = rand(3, 3);
output = L5.evaluate(input);

expected = input * scale + offset;
assert(max(abs(output(:) - expected(:))) < 1e-10, 'Scale+Offset should work correctly');

%% Test 6: ElementwiseAffineLayer with ImageStar reachability
rng(42);

scale = 2.0;
offset = 0.5;
L6 = ElementwiseAffineLayer('affine_reach', scale, offset, true, true);

% Create ImageStar input
V = zeros(3, 3, 1, 2);
V(:,:,1,1) = rand(3, 3) * 0.5;
V(:,:,1,2) = rand(3, 3) * 0.1;

input_is = ImageStar(V, [1; -1], [1; 1], -1, 1);

% Compute reachable set - just verify it produces output
output_is = L6.reach_star_single_input(input_is);
assert(~isempty(output_is), 'Should produce output ImageStar');

