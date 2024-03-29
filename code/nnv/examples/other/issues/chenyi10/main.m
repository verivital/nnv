%% Reproduce and fix error on issue (#215
% https://github.com/verivital/nnv/issues/215

net = importNetworkFromTensorFlow('my_model');
dataArray1 = rand(8,10);
dlX1 = dlarray(dataArray1,"TC");
net = initialize(net, dlX1);
summary(net)
F = matlab2nnv(net);

%% Evaluate

% Random input
IM = rand(8,10);
y = F.evaluateSequence(IM);

%% Reach

% create input set
IM = rand(8,10);
LB = IM-0.01;
UB = IM+0.01;
K = ImageStar(IM, LB, UB);

% compute reachability
reachMethod = 'approx-star'; % 使用近似星形集方法
F.reachMethod = reachMethod; % 设置网络的可达集计算方法
R = F.reachSequence(K);

