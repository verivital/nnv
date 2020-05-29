%to run this as a test, use results_cnn_image=runtests('test_cnn_image')
%requirements: file must start or end with test
%each test starts with two percent signs followed by the name
%shared vairables must appear before first test
%variables made by a test are not available to other tests.


%shared variables
% Load the trained model 
net = vgg16();

% See details of the architecture 
net.Layers

% Read the image to classify 
I0 = imread('peppers.png');

% Adjust size of the image 
sz = net.Layers(1).InputSize; 
I = I0(1:sz(1),1:sz(2),1:sz(3));

n = size(I);
N = n(1)*n(2)*n(3);


%___________________________________________________________________________________________________
%tests below originally taken from test_ImageInputLayer_evaluate.m


%% test 1: ImageInputLayer Evaluate





L1 = ImageInputLayer.parse(net.Layers(1));


Y1 = L1.evaluate(I);




%___________________________________________________________________________________________________
%tests below originally taken from test_ImageInputLayer_parse.m


%% test 2: ImageInputLayer parse



L1 = ImageInputLayer.parse(net.Layers(1));






%___________________________________________________________________________________________________
%tests below originally taken from test_ImageInputLayer_reach.m


%% test 3: ImageInputLayer reach

L1 = ImageInputLayer.parse(net.Layers(1));

I4 = reshape(I, [N,1]);
I5 = I4;
del = 200;

for i=1:N
    if I4(i) >= 255 - del
        I5(i) = 255 - I4(i);
    else
        I5(i) = 0;
    end
end

center = reshape(I4, [n(1), n(2), n(3)]); % center image matrix
basis_mat = reshape(I5, [n(1), n(2), n(3)]); % basis image matrix

C = [1;-1];   % 0% <= alpha <= bv percentage of brightening attack
bv = 0.0000005;
d = [bv; 0];
pred_lb = 0;
pred_ub = bv;

% normalized ImageStar Input Set using InputLayer of the VGG16
% note***: the InputImageLayer of VGG16 substract an image with the mean of
% images in the training set. This information is hidden, we do not know.
% Therefore, we need to use the InputImageLayer of VGG16 to compute the
% normalized ImageStar Input Set

V(:,:,:,1) = activations(net, center, 1);
V(:,:,:,2) = cast(basis_mat, 'single');
V = double(V); % use double precison for analysis

IS = ImageStar(V, C, d, pred_lb, pred_ub);

IS1 = L1.reach(IS);



%___________________________________________________________________________________________________
%tests below originally taken from test_ImageInputLayer_reach_zono.m


%% test 4: ImageInputLayer reach zono





I4 = reshape(I, [N,1]);
I5 = I4;
del = 200;

for i=1:N
    if I4(i) >= 255 - del
        I5(i) = 255 - I4(i);
    else
        I5(i) = 0;
    end
end

center = reshape(I4, [n(1), n(2), n(3)]); % center image matrix
basis_mat = reshape(I5, [n(1), n(2), n(3)]); % basis image matrix

V = double(cat(4, center, basis_mat));


tic;
image = ImageZono(V);
toc;



































