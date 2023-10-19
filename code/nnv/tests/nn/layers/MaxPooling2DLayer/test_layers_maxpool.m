%to run this as a test, use results_layers_maxpool=runtests('test_layers_maxpool')
%requirements: file must start or end with test
%each test starts with two percent signs followed by the name
%shared vairables must appear before first test
%variables made by a test are not available to other tests.


%shared variables


%___________________________________________________________________________________________________
%tests below originally taken from test_MaxPooling2DLayer_compute_maxMap.m


%% test 1: maxPooling2DLayer max Map

% this example is from (see max pooling layer)
%https://en.wikipedia.org/wiki/Convolutional_neural_network

I = [1 0 2 3; 4 6 6 8; 3 1 1 0; 1 2 2 4]; % input
L = MaxPooling2DLayer([2 2], [2 2], [0 0 0 0]);
maxMap = L.compute_maxMap(I);


I_out=[6 8; 3 4];

assert(isequal(I_out, maxMap));




%___________________________________________________________________________________________________
%tests below originally taken from test_MaxPooling2DLayer_constructor.m


%% test 2: maxPooling2DLayer constructors


% construct a Average Pooling 2D Layer object

L1 = MaxPooling2DLayer('test_max_pooling_2d_layer', [2 2], [1 1], [0 0 0 0]);
L2 = MaxPooling2DLayer();
L3 = MaxPooling2DLayer([3 3], [1 1], [0 0 0 0]);




%___________________________________________________________________________________________________
%tests below originally taken from test_MaxPooling2DLayer_evaluation.m


%% test 3: maxPooling2DLayer evaluation



% original input volume: color image with 3 channels
inputVol(:, :, 1) = [0 0 2 0; 1 2 0 2 ; 0 0 2 2; 0 2 2 2; 2 2 2 1]; % channel 1 input matrix
inputVol(:, :, 2) = [1 2 2 1; 2 1 2 0; 2 2 2 0; 1 1 1 0; 1 0 2 2]; % channel 2 input matrix
inputVol(:, :, 3) = [0 0 2 2; 0 2 1 1; 0 2 0 0; 0 2 1 0; 1 2 1 0]; % channel 3 input matrix

L = MaxPooling2DLayer([2 2], [2 2], [0 0 0 0]);

y = L.evaluate(inputVol);



outputVol(:, :, 1)=[2, 2; 2, 2];
outputVol(:, :, 2)=[2, 2; 2, 2];
outputVol(:, :, 3)=[2, 2; 2, 1];

assert(isequal(y, outputVol));



%___________________________________________________________________________________________________
%tests below originally taken from test_MaxPooling2DLayer_get_zero_padding_input.m


%% test 4: maxPooling2DLayer  zero padding





% original input volume: color image with 3 channels
inputVol(:, :, 1) = [2 0 1 2; 1 0 2 2; 1 2 2 0; 1 2 0 0; 1 0 1 1]; % channel 1 input matrix
inputVol(:, :, 2) = [0 0 1 0; 0 0 2 1; 1 1 0 1; 1 1 0 2; 2 1 2 0]; % channel 2 input matrix
inputVol(:, :, 3) = [1 2 2 1; 2 0 0 2; 0 0 1 0; 1 2 0 2; 1 0 2 1]; % channel 3 input matrix


% construct input with padding operation
paddingSize = [1 1 1 1];

L = MaxPooling2DLayer();
L.set_padding(paddingSize);

I = L.get_zero_padding_input(inputVol);

for i=1:3
    display(inputVol(:,:,i));
    display(I(:,:,i));
end


assert(isequal(inputVol(:, :, :), I(2:end-1, 2:end-1, :)));

hor_zero=zeros(1, size(I, 2), size(I, 3));
ver_zero=zeros(size(I, 1), 1, size(I, 3));

assert(isequal(hor_zero, I(1, :, :)));
assert(isequal(hor_zero, I(end, :, :)));

assert(isequal(ver_zero, I(:, 1, :)));
assert(isequal(ver_zero, I(:, end, :)));


%___________________________________________________________________________________________________
%tests below originally taken from test_MaxPooling2DLayer_reach_star_approx.m


%% test 5: maxPooling2DLayer reach star approx



% test constructor for ImageStar class
IM(:,:,1) = [1 1 0 1; 0 0 1 1; 1 0 1 0; 1 1 1 1]; % center image channel 1
IM(:,:,2) = [0 1 0 0; 1 0 0 1; 0 1 1 0; 0 0 0 1]; % center image channel 2
IM(:,:,3) = [1 1 1 1; 1 1 0 1; 0 1 1 0; 1 0 1 0]; % center image channel 3

LB(:,:,1) = [-0.1 -0.2 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]; % attack on pixel (1,1) and (1,2)
LB(:,:,2) = [-0.1 -0.15 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]; 
LB(:,:,3) = LB(:,:,2);

UB(:,:,1) = [0.1 0.2 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
UB(:,:,2) = [0.1 0.15 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
UB(:,:,3) = UB(:,:,2);

image = ImageStar(IM, LB, UB);

L = MaxPooling2DLayer([2 2], [2 2], [0 0 0 0]);

output = L.reach_star_approx(image);

% profile on;
exact_output = L.reach_star_exact(image);
% profile viewer;




%___________________________________________________________________________________________________
%tests below originally taken from test_MaxPooling2DLayer_reach_star_exact.m


%% test 6: maxPooling2DLayer reach star exact


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


%___________________________________________________________________________________________________
%tests below originally taken from test_MaxPooling2DLayer_reach_zono.m


%% test 7: maxPooling2DLayer reach zono

% test constructor for ImageStar class
LB(:,:,1) = [-0.1 -0.2 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]; % attack on pixel (1,1) and (1,2)
LB(:,:,2) = [-0.1 -0.15 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]; 
LB(:,:,3) = LB(:,:,2);

UB(:,:,1) = [0.1 0.2 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
UB(:,:,2) = [0.1 0.15 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
UB(:,:,3) = UB(:,:,2);

image = ImageZono(LB, UB);

L = MaxPooling2DLayer([2 2], [2 2], [0 0 0 0]);


%profile on;
tic;
output = L.reach_zono(image);
toc;
%profile viewer;













