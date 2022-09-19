close all; 
clear; 
clc;

dataset_ = "CIFAR10";
net_ = 'CIFAR10_FNNsmall_sigmoid';
n_ = 'FNNsmall';

net_dir = sprintf('%s/nets/%s/%s.mat', dataset_,n_,net_)
image_dir = sprintf('%s/data/%s.csv', dataset_,net_)

% reachMethod = 'approx-star';
% reachMethod = 'abs-dom'; %Star AbsDom (LP)
reachMethod = 'approx-rstar-0';
% reachMethod = 'approx-rstar-2';
% reachMethod = 'approx-rstar-4';

relaxFactor = [0];
numCores = 1;
disp_opt = 0;
lp_solver = 'linprog';
sigmoidal = 1;
normalized = 1;

%% load network
load(net_dir);
nnv_net = net2nnv_net(net, lp_solver);

%% load images
csv_data = csvread(image_dir);
IM_labels = csv_data(:,1);
IM_data = csv_data(:,2:end)';

eps = [0.004, 0.006, 0.008, 0.010, 0.012, 0.014];

N = size(IM_data, 2);
K = length(relaxFactor);
M = length(eps);

r = zeros(K, M); % percentage of images that are robust
rb = cell(K, M); % detail robustness verification result
cE = cell(K, M); % detail counterexamples
vt = cell(K, M); % detail verification time
cands = cell(K,M); % counterexample
total_vt = zeros(K, M); % total verification time

for i=1:K
    for j=1:M
        fprintf('\tepsilon: %f\n',eps(j)');
        images = attack_images(IM_data, eps(j), reachMethod, normalized); 
        labels = IM_labels + 1;
        t = tic;
        [r(i,j), rb{i,j}, cE{i,j}, cands{i,j}, vt{i,j}] = nnv_net.evaluateRBN(images(1:N), labels(1:N), reachMethod, numCores, relaxFactor , disp_opt, lp_solver, sigmoidal);
        total_vt(i,j) = toc(t);
    end
end

T = table;
ep = [];
VT = [];
RB = [];
US = [];
UK = [];
for i=1:K
    ep = [ep; eps'];
    unsafe = zeros(M,1);
    robust = zeros(M,1);
    unknown = zeros(M,1);
    for j=1:M
        unsafe(j) = sum(rb{i,j}==0);
        robust(j) = sum(rb{i,j} == 1);
        unknown(j) = sum(rb{i,j}==2);
    end
    RB = [RB; robust];
    US = [US; unsafe];
    UK = [UK; unknown];
    VT = [VT; total_vt(i,:)'];
end
T.epsilon = ep;
T.robustness = RB;
T.unsafe = US;
T.unknown = UK;
T.verifyTime = VT;

fprintf('%s', reachMethod);
T

save_ = sprintf('result/%s_%s', net_, reachMethod)
save(save_, 'lp_solver', 'T', 'r', 'rb', 'cE', 'cands', 'vt', 'total_vt');


function images = attack_images(in_images, epsilon, reachMethod, normalized)
    if normalized
        max_px = 1.0;
    else
        max_px = 255.0;
    end
    
    
    N = size(in_images, 2);
    for n = 1:N
        image = in_images(:, n);
        if normalized
            image = image/255.0;
        end
        lb = image - epsilon;
        ub = image + epsilon;
        ub(ub > max_px) = max_px;
        lb(lb < 0.0) = 0.0;
        
        if strcmp(reachMethod,'approx-star') || strcmp(reachMethod, 'abs-dom')
            images(n) = Star(lb, ub);
        elseif strcmp(reachMethod,'approx-rstar-0')
            images(n) = RStar0(lb, ub, inf);
        elseif strcmp(reachMethod,'approx-rstar-2') || strcmp(reachMethod,'approx-rstar-4')
            images(n) = RStar(lb, ub, inf);
        elseif strcmp(reachMethod,'approx-zono')
            B = Box(lb, ub);
            images(n) = B.toZono;
        else
           error('unsupported reachMethod for evaluateRBN')
        end
    end
end

function nnv_net = net2nnv_net(net, lp_solver)  
    if strcmp(net.Layers(4).Type,'Sigmoid')
    act_fn = 'logsig';
    elseif strcmp(net.Layers(4).Type,'Tanh')
        act_fn = 'tansig';
    end

    L = [];
    for i = 3:2:length(net.Layers)-4
        L1 = LayerS(net.Layers(i).Weights, net.Layers(i).Bias, act_fn);
        L1.lp_solver = lp_solver;
        L = [L L1]; 
    end
    L2 = LayerS(net.Layers(i+2).Weights, net.Layers(i+2).Bias, 'purelin');
    L2.lp_solver = lp_solver;
    nnv_net = FFNNS([L L2]);
    nnv_net.lp_solver = lp_solver;
end
