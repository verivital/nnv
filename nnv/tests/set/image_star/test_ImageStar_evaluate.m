load one_image.mat;
% constructing an imageStar using Star2D
% constructing an imageStar using Star2D
V(:,:,1) = one_image; % center matrix
V(:,:,2) = rand(28,28);  % basic matrix

% constraint: -1<= a <= 1
Constr_mat = [1; -1];
Constr_vec = [1; 1]; 
pred_lb = -1;
pred_ub = 1; 

image = ImageStar2(V, Constr_mat, Constr_vec, pred_lb, pred_ub);

image1 = image.evaluate(-1);

figure;
subplot(1,2,1);
imshow(image.V(:,:,1)); % center image

subplot(1,2,2)
imshow(image1); % sampled image from imagestar


