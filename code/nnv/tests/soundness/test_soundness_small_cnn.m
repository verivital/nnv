% test_soundness_small_cnn
% End-to-end soundness test for small CNN components
% Tests layer pipelines with soundness verification
% To run: results = runtests('test_soundness_small_cnn')

%% Test 1: Conv -> ReLU pipeline soundness
% Test a simple Conv2D -> ReLU pipeline
rng(42);

% Layer 1: Conv2D (3x3 filter, 1->1 channel)
W1 = randn(3, 3, 1, 1) * 0.5;
b1 = 0.1;
L1 = Conv2DLayer(W1, b1);

% Layer 2: ReLU
L2 = ReluLayer();

% Create input ImageStar (6x6x1 with 2 predicate variables)
V = zeros(6, 6, 1, 3);
V(:,:,1,1) = rand(6, 6);  % center
V(2,2,1,2) = 0.1;         % basis 1
V(4,4,1,3) = 0.1;         % basis 2

C = [eye(2); -eye(2)];
d = [1; 1; 1; 1];
pred_lb = [-1; -1];
pred_ub = [1; 1];

input_is = ImageStar(V, C, d, pred_lb, pred_ub);

% Propagate through Conv -> ReLU
IS1 = L1.reach(input_is, 'exact-star');
IS2 = L2.reach(IS1, 'approx-star');

% Convert to cell array for iteration
if ~iscell(IS2)
    if length(IS2) > 1
        output_sets = cell(1, length(IS2));
        for idx = 1:length(IS2)
            output_sets{idx} = IS2(idx);
        end
    else
        output_sets = {IS2};
    end
else
    output_sets = IS2;
end

% Verify soundness: sample inputs and check outputs
n_samples = 30;
for i = 1:n_samples
    alpha = pred_lb + (pred_ub - pred_lb) .* rand(2, 1);
    if any(C * alpha > d)
        continue;
    end

    input_concrete = V(:,:,:,1) + alpha(1)*V(:,:,:,2) + alpha(2)*V(:,:,:,3);
    x = L1.evaluate(input_concrete);
    output_concrete = L2.evaluate(x);

    contained = false;
    for j = 1:length(output_sets)
        if soundness_test_utils.verify_imagestar_containment(output_sets{j}, output_concrete, 1e-4)
            contained = true;
            break;
        end
    end
    assert(contained, 'Conv-ReLU pipeline soundness violation at sample %d', i);
end

%% Test 2: Conv -> AvgPool pipeline soundness
rng(42);

W1b = randn(2, 2, 1, 1) * 0.5;
b1b = 0;
L1b = Conv2DLayer(W1b, b1b);

L2b = AveragePooling2DLayer([2 2], [2 2], [0 0 0 0]);

% 4x4 input
Vb = zeros(4, 4, 1, 2);
Vb(:,:,1,1) = rand(4, 4);
Vb(2,2,1,2) = 0.2;

input_isb = ImageStar(Vb, [1; -1], [0.5; 0.5], -0.5, 0.5);

% Propagate
IS1b = L1b.reach(input_isb, 'exact-star');
IS2b = L2b.reach(IS1b, 'exact-star');

% Verify soundness
for i = 1:30
    alpha = -0.5 + rand();
    input_concrete = Vb(:,:,:,1) + alpha * Vb(:,:,:,2);
    x = L1b.evaluate(input_concrete);
    output_concrete = L2b.evaluate(x);

    contained = soundness_test_utils.verify_imagestar_containment(IS2b, output_concrete, 1e-5);
    assert(contained, 'Conv-AvgPool pipeline soundness violation at sample %d', i);
end

%% Test 3: Multi-layer linear pipeline (Conv -> AvgPool -> Conv)
rng(42);

W1c = randn(2, 2, 1, 1) * 0.5;
b1c = 0.1;
L1c = Conv2DLayer(W1c, b1c);

L2c = AveragePooling2DLayer([2 2], [2 2], [0 0 0 0]);

W3c = randn(2, 2, 1, 1) * 0.5;
b3c = -0.1;
L3c = Conv2DLayer(W3c, b3c);

% 6x6 input -> Conv(2x2) -> 5x5 -> Pool(2x2,s2) -> 2x2 -> Conv(2x2) -> 1x1
Vc = zeros(6, 6, 1, 3);
Vc(:,:,1,1) = rand(6, 6);
Vc(1,1,1,2) = 0.1;
Vc(4,4,1,3) = 0.1;

Cc = [eye(2); -eye(2)];
dc = ones(4, 1);
input_isc = ImageStar(Vc, Cc, dc, [-1;-1], [1;1]);

% Propagate
IS1c = L1c.reach(input_isc, 'exact-star');
IS2c = L2c.reach(IS1c, 'exact-star');
IS3c = L3c.reach(IS2c, 'exact-star');

% Verify corners
corners = [-1 -1; -1 1; 1 -1; 1 1];
for i = 1:4
    alpha = corners(i, :)';

    input_concrete = Vc(:,:,:,1) + alpha(1)*Vc(:,:,:,2) + alpha(2)*Vc(:,:,:,3);
    x = L1c.evaluate(input_concrete);
    x = L2c.evaluate(x);
    output_concrete = L3c.evaluate(x);

    contained = soundness_test_utils.verify_imagestar_containment(IS3c, output_concrete, 1e-5);
    assert(contained, 'Linear pipeline corner %d not contained', i);
end

fprintf('All small CNN soundness tests passed!\n');
