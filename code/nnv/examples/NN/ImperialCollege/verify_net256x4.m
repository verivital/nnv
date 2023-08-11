% Load the network
load net256x4.mat;

% Load the data
load inputSet.mat;

% Select number of sets to evaluate (max of 50)
N = 12;
rob_res2 = zeros(N,2);
reachOptions = struct;
reachOptions.reachMethod = 'approx-star';

% verify the network with eps = 0.02
t = tic;
parfor i=1:N
    rob_res2(i,1) = net.verify_robustness(S_eps_002(i), reachOptions, labels(i));
end
verify_time2 = toc(t);

% verify the network with eps = 0.05
rob_res5 = zeros(N,2);
t = tic;
parfor i=1:N
    rob_res5(i,1) = net.verify_robustness(S_eps_005(i), reachOptions, labels(i));
end
verify_time5 = toc(t);


% Results summary
epsilon = [0.02; 0.05];
verify_time = [verify_time2; verify_time5];
safe = [sum(rob_res2(:,1)==1); sum(rob_res5(:,1) == 1)];
unsafe = [sum(rob_res2(:,1) == 0); sum(rob_res5(:,1) == 0)];
unknown = [sum(rob_res2(:,1) == 2); sum(rob_res5(:,1) == 2)];

T = table(epsilon, safe, unsafe, unknown, verify_time)

