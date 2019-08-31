load one_image.mat
% constructing an imageStar using Star2D
Center = one_image; % center matrix
Basis = rand(28,28); % basic matrix

V(:,:,1,1) = Center;
V(:,:,1,2) = Basis;

% constraint: -1<= a <= 1
Constr_mat = [1; -1];
Constr_vec = [1; 1]; 
pred_lb = -1;
pred_ub = 1; 

input_image = ImageStar(V, Constr_mat, Constr_vec, pred_lb, pred_ub);

L = AveragePooling2DLayer([6 4], [4 4], [1 1 0 0]);

output_image = L.reach_star_single_input(input_image);

sampled_images = output_image.sample(2);


figure;
subplot(1,3,1);
imshow(input_image.V(:,:,1)); % center image
title('28x28 input image');
subplot(1,3,2);
imshow(sampled_images{1, 1}); % the first sampled image from imagestar
title('7x7 1st output image');
subplot(1,3,3);
imshow(sampled_images{1, 2}); % the second sampled image from imagestar
title('7x7 2nd output image');


