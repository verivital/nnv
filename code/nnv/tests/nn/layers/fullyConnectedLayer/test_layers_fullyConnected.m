%to run this as a test, use results_layers_fullyConnected=runtests('test_layers_fullyConnected')
%requirements: file must start or end with test
%each test starts with two percent signs followed by the name
%shared vairables must appear before first test
%variables made by a test are not available to other tests.


%shared variables
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

image_star = ImageStar(IM, LB, UB);
image_zono = ImageZono(LB, UB);


%___________________________________________________________________________________________________
%tests below originally taken from test_FullyConnectedLayer_constructor.m


%% test 1: FullyConnected Constructor


% construct FullyConnectedLayer objects
W = rand(4, 48);
b = rand(4,1);
fc = FullyConnectedLayer(W,b);
fc1 = FullyConnectedLayer('fc1', W, b);



%___________________________________________________________________________________________________
%tests below originally taken from test_FullyConnectedLayer_reach_star_exact.m


%% test 2: FullyConnected reach star exact



% construct a FullyConnectedLayer object
W = rand(4, 48);
b = rand(4,1);
fc = FullyConnectedLayer(W,b);

%output = fc.reach_star_exact(image_star);
%output2 = fc.reach_zono_exact(image_zono);

output = fc.reach_star_single_input(image_star);
output2 = fc.reach_zono(image_zono);

%___________________________________________________________________________________________________
%tests below originally taken from test_FullyConnectedLayer_reach_zono_exact.m


%% test 3: FullyConnected reach zono exact



% construct a FullyConnectedLayer object
W = rand(4, 48);
b = rand(4,1);
fc = FullyConnectedLayer(W,b);

%output = fc.reach_zono_exact(image_zono);
output = fc.reach_zono(image_zono);








































