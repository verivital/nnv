% Load data and network
load('../../../io/models/triangle_net.mat');
nnvSegNet = matlab2nnv(net);

% Load data
dataSetDir = fullfile(toolboxdir('vision'),'visiondata','triangleImages');
imageDir = fullfile(dataSetDir,'trainingImages');
labelDir = fullfile(dataSetDir,'trainingLabels');
imds = imageDatastore(imageDir);

I0 = readimage(imds, 1);

% Adjust size of the image 
sz = net.Layers(1).InputSize; 
I = single(I0(1:sz(1),1:sz(2),1:sz(3)));

n = size(I);
N = numel(I);

L1 = ImageInputLayer.parse(net.Layers(1));

%% test 1: ImageInputLayer Evaluate
Y1 = L1.evaluate(I);

%% test 3: ImageInputLayer reach

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

center = reshape(I4, n); % center image matrix
basis_mat = reshape(I5, n); % basis image matrix

C = [1;-1];   % 0% <= alpha <= bv percentage of brightening attack
bv = 0.0000005;
d = [bv; 0];
pred_lb = 0;
pred_ub = bv;

V(:,:,:,1) = activations(net, center, net.Layers(1).Name);
V(:,:,:,2) = cast(basis_mat, 'single');
V = double(V); % use double precison for analysis

IS = ImageStar(V, C, d, pred_lb, pred_ub);
IS1 = L1.reach(IS, 'approx-star');

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

center = reshape(I4, n); % center image matrix
basis_mat = reshape(I5, n); % basis image matrix

V = double(cat(4, center, basis_mat));
image = ImageZono(V);

IS2 = L1.reach(image, 'approx-zono');