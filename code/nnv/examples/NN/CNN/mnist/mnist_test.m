%I0 = imread('peppers.png');

opt_retrain = 0;

if opt_retrain
    o_dir = pwd;
    openExample('nnet/TrainABasicConvolutionalNeuralNetworkForClassificationExample');
    TrainABasicConvolutionalNeuralNetworkForClassificationExample
    cd(o_dir)
end

nnvNet = CNN.parse(net, 'mnist_builtin');

i_test = ceil(rand(1,1) * length(imds.Files));
figure ;
image_test_path = imds.Files{i_test} ; 
image_test = imread(image_test_path) ; 
%image_test = imrotate(image_test,180) ; 
imshow(image_test) ; 
v_test = net.predict(image_test); 
[i_test_v,i_test_c] = max(v_test) ; 
test_label = net.Layers(end).Classes(i_test_c) ; 
real_label = imds.Labels(i_test) ; 
{test_label == real_label,test_label,real_label} % show if classification prediction is correct or not

I0 = image_test;


% close all ; i_test = ceil(rand(1,1) * length(imds.Files)); figure ; image_test_path = imds.Files{i_test} ; image_test = imread(image_test_path) ; image_test = imrotate(image_test,180) ; imshow(image_test) ; v_test = net.predict(image_test); [i_test_v,i_test_c] = max(v_test) ; test_label = net.Layers(end).Classes(i_test_c) ; real_label = imds.Labels(i_test) ; {test_label == real_label,test_label,real_label}

% Adjust size of the image 
sz = nnvNet.Layers{1}.InputSize; 
I = I0(1:sz(1),1:sz(2),1:sz(3));

n = size(I);
N = n(1)*n(2);

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

center = reshape(I4, [n(1), n(2)]); % center image matrix
basis_mat = reshape(I5, [n(1), n(2)]); % basis image matrix

C = [1;-1];   % 0% <= alpha <= bv percentage of brightening attack
bv = 0.0000001;
d = [bv; 0];
pred_lb = 0;
pred_ub = bv;

V(:,:,:,1) = center;
V(:,:,:,2) = basis_mat;

IS = ImageStar(V, C, d, pred_lb, pred_ub);

fprintf('\n========= PARSE VGG16 FOR REACHABILITY ANALYSIS ============\n');

fprintf('\n======= DO REACHABILITY ANLAYSIS WITH EXACT-STAR METHOD ======\n');

numCores = 1;
nnvNet.reach(IS, 'exact-star', numCores);


figure; hold on;
plot([0:9], squeeze(nnvNet.reachSet{end}.V(:,:,:,1)));
plot([0:9], squeeze(nnvNet.reachSet{end}.V(:,:,:,2)));

