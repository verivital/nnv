load one_image.mat
% constructing an imageStar using Star2D
Center = one_image; % center matrix
Basis = rand(28,28); % basic matrix

V(:,:,1) = Center;
V(:,:,2) = Basis;

% constraint: -1<= a <= 1
Constr_mat = [1; -1];
Constr_vec = [1; 1]; 
pred_lb = -1;
pred_ub = 1; 

image = ImageStar(V, Constr_mat, Constr_vec, pred_lb, pred_ub);

image1 = image.sample(2);

figure;
subplot(1,3,1);
imshow(image.V(:,:,1)); % center image
subplot(1,3,2);
imshow(image1{1,1}); % the first sampled image from imagestar
subplot(1,3,3);
imshow(image1{1,2}); % the second sampled image from imagestar


