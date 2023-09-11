fprintf('\n\n=============================LOAD VGG19 ======================\n');

% Load the trained model 
net = vgg19();

fprintf('\n\n======================== PARSING VGG19 =======================\n');
nnvNet = matlab2nnv(net);


fprintf('\n\n=========CONSTRUCT INPUT SET (AN IMAGESTAR SET) =============\n');
load image_data.mat;
V(:,:,:,1) = double(ori_image);
V(:,:,:,2) = double(dif_image);

pred_lb = 0.5;
delta = 0.0000001;
pred_ub = pred_lb + delta;
C = [1;-1];   % pred_lb % <= alpha <= pred_ub percentage of FGSM attack
d = [pred_ub; -pred_lb];
IS = ImageStar(double(V), C, d, pred_lb, pred_ub);

fprintf('\n\n======= DO REACHABILITY ANLAYSIS WITH EXACT-STAR METHOD ======\n');

reachOptions = struct;
reachOptions.reachMethod = 'exact-star';
nnvNet.reach(IS, reachOptions);

exactReachSet = nnvNet.reachSet{end};

fprintf('\n\n========REACHABILITY IS DONE IN %.5f SECONDS==========\n', sum(nnvNet.reachTime));

% Get the reachTime per layer
RT = nnvNet.reachTime;
conv_rT = 0;
fc_rT = 0;
maxpool_rT = 0;
relu_rT = 0;
for i = 1:length(nnvNet.Layers)
    % Convolutional layers
    if isa(nnvNet.Layers{i}, 'Conv2DLayer')
        conv_rT = conv_rT + nnvNet.reachTime(i);
    elseif isa(nnvNet.Layers{i}, 'FullyConnectedLayer')
        fc_rT = fc_rT + nnvNet.reachTime(i);
    elseif isa(nnvNet.Layers{i}, 'MaxPooling2DLayer')
        maxpool_rT = maxpool_rT + nnvNet.reachTime(i);
    elseif isa(nnvNet.Layers{i}, 'ReluLayer')
        relu_rT = relu_rT + nnvNet.reachTime(i);
    end
end

% Visualize results
figure;
c = categorical({'Convolutional 2D (16)','Fully Connected (3)', 'Max Pooling (5)', 'ReLU (18)'});
reachTime = [conv_rT fc_rT maxpool_rT relu_rT];
bar(c,reachTime);
text(1:length(reachTime),reachTime,num2str(reachTime'),'vert','bottom','horiz','center');

