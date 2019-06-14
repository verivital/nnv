load one_image.mat
% constructing an imageStar using Star2D
Center = one_image; % center matrix
Basis = rand(28,28); % basic matrix

V = cell(1, 2);
V{1} = Center;
V{2} = Basis;

% constraint: -1<= a <= 1
Constr_mat = [1; -1];
Constr_vec = [1; 1]; 
pred_lb = -1;
pred_ub = 1; 

image = ImageStar(V, Constr_mat, Constr_vec, pred_lb, pred_ub);

image1 = image.evaluate(-1);

figure; 
imshow(image.V{1}{1}); % center image
%title('Center image of the imagestar');

figure;
imshow(image1); % sampled image from imagestar
%title('Sampled image from the imagestar');

