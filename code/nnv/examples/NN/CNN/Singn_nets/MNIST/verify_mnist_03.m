
%% load network and input set
load mnist03.mat;

%% verify the image

numCores = 6; 
N = 100; % number of Images tested maximum is 100
eps = [0.3];
relaxFactor = [0]; 
% eps = [0.1 0.15 0.2 0.25 0.3]; % epsilon for the disturbance
% % timeout for eps = 0.25
% relaxFactor = [0 0.25 0.5 0.75 1]; % relax factor
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
        [IS, Labels] = getInputSet(eps(j)); 
        t = tic;
        [r(i,j), rb{i,j}, cE{i,j}, cands{i,j}, vt{i,j}] = net.evaluateRBN(IS(1:N), Labels(1:N), 'approx-star', numCores, relaxFactor(i));
        total_vt(i,j) = toc(t);
    end
end

%% print the result
T = table;
rf = [];
ep = [];
VT = [];
RB = [];
US = [];
UK = [];
for i=1:K
    rf = [rf; relaxFactor(i)*ones(M,1)];
    ep = [ep; eps'];
    unsafe = zeros(M,1);
    robust = zeros(M,1);
    unknown = zeros(M,1);
    for j=1:M
        unsafe(i) = sum(rb{i,j}==0);
        robust(i) = sum(rb{i,j} == 1);
        unknown(i) = sum(rb{i,j}==2);
    end
    RB = [RB; robust];
    US = [US; unsafe];
    UK = [UK; unknown];
    VT = [VT; total_vt(i,:)'];
end
T.relaxFactor = rf; 
T.epsilon = ep;
T.robustness = RB;
T.unsafe = US;
T.unknown = UK;
T.verifyTime = VT;
T

save mnist03_results.mat T r rb cE cands vt total_vt;