% Load the network
load net256x6.mat;

% Load the data
load inputSet.mat;

% Select number of sets to evaluate (max of 50)
N = 12;
rob_res = zeros(N,2);
reachOptions = struct;
reachOptions.reachMethod = 'approx-star';

% verify the network with eps = 0.02
t = tic;
parfor i=1:N
    rob_res(i,1) = net.verify_robustness(S_eps_002(i), reachOptions, labels(i));
end
verify_time = toc(t);

% Results summary
epsilon = 0.02;
safe = sum(rob_res(:,1) == 1);
unsafe = sum(rob_res(:,1) == 0);
unknown = sum(rob_res(:,1) == 2);

T = table(epsilon, safe, unsafe, unknown, verify_time)
