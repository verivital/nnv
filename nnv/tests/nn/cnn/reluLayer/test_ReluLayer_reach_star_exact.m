% construct a FullyConnectedLayer object
L = ReluLayer();
V(:, :, 1, 1) = [-1 1; 0 2]; % channel 1 input matrix
basis = zeros(2,2);
basis(1, 1) = 1;
basis(2, 1) = 1;
V(:, :, 1, 2) = basis;

C = [1; -1];
d = [2; 2];
pred_lb = -2;
pred_ub = 2;

in_image = ImageStar(V, C, d, pred_lb, pred_ub);

images = L.reach(in_image, 'exact-star');