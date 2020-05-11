V(:, :, 1, 1) = [0 4 1 2; 2 3 2 3 ; 1  3 1 2; 2 1 3 2]; % channel 1 input matrix
basis = zeros(4,4);
basis(1, 2) = 1;
basis(4,1) = 1;
V(:, :, 1, 2) = basis;

C = [1; -1];
d = [2; 2];
pred_lb = -2;
pred_ub = 2;

in_image = ImageStar(V, C, d, pred_lb, pred_ub);

L = MaxPooling2DLayer([2 2], [2 2], [0 0 0 0]);
L1 = MaxUnpooling2DLayer;
L1.PairedMaxPoolingName = L.Name;

IS1 = L.reach(in_image, 'exact-star');
OS1 = L1.reach(IS1, 'exact-star');

IS2 = L.reach(in_image, 'approx-star');
OS2 = L1.reach(IS2, 'approx-star');