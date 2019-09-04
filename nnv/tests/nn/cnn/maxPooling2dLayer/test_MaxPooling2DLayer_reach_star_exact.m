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

images = L.reach_star_exact(in_image);

fprintf('\nFirst ImageStar:\n');
display(images(1).V)
display(images(1).C)
display(images(1).d)

fprintf('\nSecond ImageStar:\n');
display(images(2).V)
display(images(2).C)
display(images(2).d)

fprintf('\nThird ImageStar:\n');
display(images(1).V)
display(images(1).C)
display(images(1).d)
