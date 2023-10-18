%to run this as a test, use results_layers_relu=runtests('test_layers_relu')
%requirements: file must start or end with test
%each test starts with two percent signs followed by the name
%shared vairables must appear before first test
%variables made by a test are not available to other tests.


%shared variables


%___________________________________________________________________________________________________
%tests below originally taken from test_ReluLayer_constructor.m


%% test 1: fully connected layer
% construct FullyConnectedLayer objects
L = ReluLayer();
L1 = ReluLayer('relu1');




%___________________________________________________________________________________________________
%tests below originally taken from test_ReluLayer_evaluate.m


%% test 2: relu evaluate

% construct a FullyConnectedLayer object
rl = ReluLayer();
% image input set
IM(:,:,1) = [-1 1 0 -1; 0 0 1 -1; 1 0 -1 0; 1 -1 -1 1]; % center image channel 1
IM(:,:,2) = [0 1 0 0; 1 0 0 -1; 0 1 -1 0; 0 1 0 -1]; % center image channel 2
IM(:,:,3) = [1 -1 1 1; 1 -1 0 1; 0 1 -1 0; 1 0 -1 0]; % center image channel 3

output = rl.evaluate(IM);


IM_out(:,:,1) = [0 1 0 0; 0 0 1 0; 1 0 0 0; 1 0 0 1]; % center image channel 1
IM_out(:,:,2) = [0 1 0 0; 1 0 0 0; 0 1 0 0; 0 1 0 0]; % center image channel 2
IM_out(:,:,3) = [1 0 1 1; 1 0 0 1; 0 1 0 0; 1 0 0 0]; % center image channel 3

assert(isequal(output, IM_out));




%___________________________________________________________________________________________________
%tests below originally taken from test_ReluLayer_reach_star_approx.m


%% test 3: relu reach star approx


% construct a FullyConnectedLayer object
L = ReluLayer();
% image input set
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
tic;
%images = L.reach(in_image, 'approx-star');
images = L.reach(image, 'approx-star');%above line replaced....maybe it was a typo?
toc;



%___________________________________________________________________________________________________
%tests below originally taken from test_ReluLayer_reach_star_exact.m


%% test 4: relu reach star exact



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







%___________________________________________________________________________________________________
%tests below originally taken from test_ReluLayer_reach_zono.m


%% test 5: relu reach zono



% construct a FullyConnectedLayer object
L = ReluLayer();
LB(:,:,1) = [-0.1 -0.2 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]; % attack on pixel (1,1) and (1,2)
LB(:,:,2) = [-0.1 -0.15 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]; 
LB(:,:,3) = LB(:,:,2);

UB(:,:,1) = [0 0.2 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
UB(:,:,2) = [0.1 0.15 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
UB(:,:,3) = UB(:,:,2);

image = ImageZono(LB, UB);

%image1 = L.reach_zono(image, 'single');
image1 = L.reach_zono(image);














