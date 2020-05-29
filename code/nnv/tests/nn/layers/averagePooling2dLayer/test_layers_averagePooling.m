%to run this as a test, use results_cnn_averagePooling=runtests('test_cnn_averagePooling')
%requirements: file must start or end with test
%each test starts with two percent signs followed by the name
%shared vairables must appear before first test
%variables made by a test are not available to other tests.


%shared variables




% this example is from (see max pooling layer)
%https://en.wikipedia.org/wiki/Convolutional_neural_network


%___________________________________________________________________________________________________
%tests below originally taken from test_AveragePooling2DLayer_compute_averageMap.m

%% test 1: AveragePooling2DLayer Compute averageMap
I = [1 0 2 3; 4 6 6 8; 3 1 1 0; 1 2 2 4]; % input
L = AveragePooling2DLayer([2 2], [2 2], [0 0 0 0]);
averageMap = L.compute_averageMap(I);
display(averageMap);
display(I);

checker=[(I(1, 1)+I(1, 2)+I(2, 1)+I(2, 2))/4, (I(1, 3)+I(1, 4)+I(2, 3)+I(2, 4))/4; (I(3, 1)+I(3, 2)+I(4, 1)+I(4, 2))/4, (I(3, 3)+I(3, 4)+I(4, 3)+I(4, 4))/4];
assert(isequal(checker, averageMap))



%___________________________________________________________________________________________________
%tests originally taken from test_AveragePooling2DLayer_constructor.m


%% test 2: AveragePooling2DLayer constructor
L1 = AveragePooling2DLayer('test_average_pooling_2d_layer', [2 2], [1 1], [0 0 0 0]);
L2 = AveragePooling2DLayer();
L3 = AveragePooling2DLayer([3 3], [1 1], [0 0 0 0]);



%___________________________________________________________________________________________________
%tests below originally taken from test_AveragePooling2DLayer_evaluation.m

%% test 3: AveragePooling2DLayer evaluation
% original input volume: color image with 3 channels
inputVol(:, :, 1) = [0 0 2 0 0; 1 2 0 2 0; 0 0 2 2 0; 0 2 2 2 2; 2 2 2 1 1]; % channel 1 input matrix
inputVol(:, :, 2) = [1 2 2 1 2; 2 1 2 0 2; 2 2 2 0 1; 1 1 1 0 0; 1 0 2 2 1]; % channel 2 input matrix
inputVol(:, :, 3) = [0 0 2 2 1; 0 2 1 1 2; 0 2 0 0 1; 0 2 1 0 1; 1 2 1 0 0]; % channel 3 input matrix

L = AveragePooling2DLayer([3 3], [2 2], [0 0 0 0]);
y = L.evaluate(inputVol);

for i=1:3
    display(y(:,:,i));
    display(inputVol(:,:,i));
end




%___________________________________________________________________________________________________
%tests below originally taken from test_AveragePooling2DLayer_get_zero_padding_input.m

%% test 4: AveragePooling2DLayer get zero padding input

% original input volume: color image with 3 channels
inputVol(:, :, 1) = [2 0 1 2 1; 1 0 2 2 2; 1 2 2 0 2; 1 2 0 0 1; 1 0 1 1 2]; % channel 1 input matrix
inputVol(:, :, 2) = [0 0 1 0 1; 0 0 2 1 1; 1 1 0 1 1; 1 1 0 2 2; 2 1 2 0 0]; % channel 2 input matrix
inputVol(:, :, 3) = [1 2 2 1 0; 2 0 0 2 0; 0 0 1 0 1; 1 2 0 2 0; 1 0 2 1 0]; % channel 3 input matrix


% construct input with padding operation
paddingSize = [1 1 1 1];

L = AveragePooling2DLayer();
L.set_padding(paddingSize);


I = L.get_zero_padding_input(inputVol);

assert(isequal(I(2:end-1, 2:end-1, :), inputVol))

corners=[I(1, 1, :), I(1, end, :); I(end, 1, :), I(end, end, :)];
top_bot=[I(1, 2:end-1, :); I(end, 2:end-1, :)];
left_right=[I(2:end-1, 1, :), I(2:end-1, end, :)];

assert(isempty(find(corners)))
assert(isempty(find(top_bot)))
assert(isempty(find(left_right)))


%___________________________________________________________________________________________________
%tests below originally taken from test_AveragePooling2DLayer_reach_star_digit_one_example.m


%% test 5: AveragePooling2DLayer reach star digit one example
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


%TODO: figure out if there is something we need to check the output.



%___________________________________________________________________________________________________
%tests below originally taken from test_AveragePooling2DLayer_reach_star_exact_single_input.m


%% test 6: AveragePooling2DLayer reach star exact single input
% original input volume: color image with 3 channels
IM(:, :, 1) = [0 0 2 0 0; 1 2 0 2 0; 0 0 2 2 0; 0 2 2 2 2; 2 2 2 1 1]; % channel 1 input matrix
IM(:, :, 2) = [1 2 2 1 2; 2 1 2 0 2; 2 2 2 0 1; 1 1 1 0 0; 1 0 2 2 1]; % channel 2 input matrix
IM(:, :, 3) = [0 0 2 2 1; 0 2 1 1 2; 0 2 0 0 1; 0 2 1 0 1; 1 2 1 0 0]; % channel 3 input matrix


LB(:,:,1) = [-0.1 -0.2 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0]; % attack on pixel (1,1) and (1,2)
LB(:,:,2) = [-0.1 -0.15 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0]; 
LB(:,:,3) = LB(:,:,2);

UB(:,:,1) = [0.1 0.2 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0;0 0 0 0 0];
UB(:,:,2) = [0.1 0.15 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0;0 0 0 0 0];
UB(:,:,3) = UB(:,:,2);

input = ImageStar(IM, LB, UB);

L = AveragePooling2DLayer([3 3], [2 2], [0 0 0 0]);

Y = L.reach_star_single_input(input);

%todo, figure out what to do for checking this

