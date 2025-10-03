% Test SoftmaxLayer functionality
% To run: results = runtests('test_SoftmaxLayer')

%% Test 1: SoftmaxLayer constructor - default
L = SoftmaxLayer();
assert(strcmp(L.Name, 'SoftmaxLayer'));
assert(L.NumInputs == 1);
assert(L.NumOutputs == 1);

%% Test 2: SoftmaxLayer constructor - with name
L = SoftmaxLayer('softmax1');
assert(strcmp(L.Name, 'softmax1'));

%% Test 3: SoftmaxLayer constructor - with all parameters
L = SoftmaxLayer('softmax2', 1, 1, {'input'}, {'output'});
assert(strcmp(L.Name, 'softmax2'));
assert(strcmp(L.InputNames{1}, 'input'));
assert(strcmp(L.OutputNames{1}, 'output'));

%% Test 4: SoftmaxLayer evaluate - vector input
L = SoftmaxLayer();

% Test with simple vector
input = [1; 2; 3; 4; 5];
output = L.evaluate(input);

% Check softmax properties: sum = 1, all positive
assert(abs(sum(output) - 1) < 1e-6, 'Softmax outputs should sum to 1');
assert(all(output > 0), 'All softmax outputs should be positive');

% Check that higher inputs produce higher outputs
assert(output(5) > output(1), 'Softmax should preserve order');

%% Test 5: SoftmaxLayer evaluate - 3D image
L = SoftmaxLayer();

% Create 3D image (H x W x C)
IM(:,:,1) = [1 2; 3 4];
IM(:,:,2) = [2 3; 4 5];
IM(:,:,3) = [3 4; 5 6];

output = L.evaluate(IM);

% Output should have same size as input
assert(all(size(output) == size(IM)), 'Output size should match input');

% Each spatial location should sum to 1 across channels
for i = 1:size(IM, 1)
    for j = 1:size(IM, 2)
        channel_sum = sum(output(i, j, :));
        assert(abs(channel_sum - 1) < 1e-6, 'Softmax should sum to 1 across channels');
    end
end

%% Test 6: SoftmaxLayer evaluate - numerical stability
L = SoftmaxLayer();

% Test with large values (numerical stability check)
input = [100; 101; 102];
output = L.evaluate(input);

% Should still sum to 1 and be all positive
assert(abs(sum(output) - 1) < 1e-6, 'Softmax should handle large values');
assert(all(output > 0), 'All outputs should be positive');

%% Test 7: SoftmaxLayer reach with ImageStar
L = SoftmaxLayer();

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
output_star = L.reach(image_star);

% Softmax layer passes through unchanged in reachability
assert(isa(output_star, 'ImageStar'), 'reach should return ImageStar');
assert(isequal(output_star.V, image_star.V), 'Softmax reach should return input unchanged');

%% Test 8: SoftmaxLayer reach with different methods
L = SoftmaxLayer();

% Create simple ImageStar
IM = rand(2, 2, 3);
LB = -0.1 * ones(2, 2, 3);
UB = 0.1 * ones(2, 2, 3);
image_star = ImageStar(IM, LB, UB);

% Test different reachability methods (all should pass through)
output1 = L.reach(image_star, 'exact-star');
output2 = L.reach(image_star, 'approx-star');

assert(isequal(output1.V, image_star.V), 'exact-star should pass through');
assert(isequal(output2.V, image_star.V), 'approx-star should pass through');

%% Test 9: SoftmaxLayer toGPU
L = SoftmaxLayer();
L_gpu = L.toGPU();
assert(isa(L_gpu, 'SoftmaxLayer'));

%% Test 10: SoftmaxLayer changeParamsPrecision
L = SoftmaxLayer();
L_single = L.changeParamsPrecision('single');
assert(isa(L_single, 'SoftmaxLayer'));
