% Load the network
load net256x2.mat;

% Load the data
load inputSet.mat;

%% Part 1. Verification (approx analysis)

% Select number of sets to evaluate (max of 50)
N = 25;
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
verify_time_approx = [verify_time2; verify_time5];
safe = [sum(rob_res2(:,1)==1); sum(rob_res5(:,1) == 1)];
unsafe = [sum(rob_res2(:,1) == 0); sum(rob_res5(:,1) == 0)];
unknown = [sum(rob_res2(:,1) == 2); sum(rob_res5(:,1) == 2)];

disp("Summary of results after step 1");

T = table(epsilon, safe, unsafe, unknown, verify_time_approx);
disp(T);

%% Part 2 - Falsification

% Results show lots of unknowns 24/25 in epsilon = 0.005)
% Let's try to find counter examples for those

% Start with those not safe with epsilon = 0.02
unk_idxs = find(rob_res2(:,1)==2)';
n_samples = 1000; % simulate the network with n_samples examples to find counter examples
t = tic;
for i=unk_idxs
    % Define unsafe region
    U = net.robustness_set(labels(i), 'max');
    % Falsify using simulations
    counter_exs = net.falsify(S_eps_002(i), U, n_samples);
    if ~isempty(counter_exs)
        rob_res2(i,1) = 0; % violated (counter example found)
    end
end
falsify_time2 = toc(t);

% Start with those not safe with epsilon = 0.02
unk_idxs = find(rob_res5(:,1)==2)';
n_samples = 200; % simulate the network with n_samples examples to find counter examples
t = tic;
for i=unk_idxs
    % Define unsafe region
    U = net.robustness_set(labels(i), 'max');
    % Falsify using simulations
    counter_exs = net.falsify(S_eps_005(i), U, n_samples);
    if ~isempty(counter_exs)
        rob_res5(i,1) = 0; % violated (counter example found)
    end
end
falsify_time5 = toc(t);

% Results summary
epsilon = [0.02; 0.05];
verify_time_approx = [verify_time2; verify_time5];
falsify_time = [falsify_time2; falsify_time5];
safe = [sum(rob_res2(:,1)==1); sum(rob_res5(:,1) == 1)];
unsafe = [sum(rob_res2(:,1) == 0); sum(rob_res5(:,1) == 0)];
unknown = [sum(rob_res2(:,1) == 2); sum(rob_res5(:,1) == 2)];

disp("Updated table after falsification attempts...");

T = table(epsilon, safe, unsafe, unknown, verify_time_approx, falsify_time);
disp(T);

%% Part 3 - Verification with exact analysis
% Results still show lots of unknown input sets
% Let's refine the verification process using the exact analysis
% (just do one of them as an example)

reachOptions = struct;
reachOptions.reachMethod = 'exact-star';
numCores = feature('numcores');
reachOptions.numCores = numCores;

% Exact-reachability 
unk_idxs = find(rob_res2(:,1)==2)';
t = tic;
rob_res2(unk_idxs(end),1) = net.verify_robustness(S_eps_002(unk_idxs(end)), reachOptions, labels(unk_idxs(end)));
toc(t);

% Could also run the evaluaiton on all  unknown images by uncommenting the
% following code:

% verify the network with eps = 0.02 (only unknowns)
% unk_idxs = find(rob_res2(:,1)==2)';
% t = tic;
% for i=unk_idxs
%     rob_res2(unk_idxs(i),1) = net.verify_robustness(S_eps_002(unk_idxs(i)), reachOptions, labels(unk_idxs(i)));
% end
% exact_time2 = toc(t);
% 
% % verify the network with eps = 0.05 (only unknowns)
% unk_idxs = find(rob_res5(:,1)==2)';
% t = tic;
% for i=unk_idxs
%     rob_res5(i,1) = net.verify_robustness(S_eps_005(i), reachOptions, labels(i));
% end
% exact_time5 = toc(t);
% 
% % Results summary
% epsilon = [0.02; 0.05];
% verify_time_approx = [verify_time2; verify_time5];
% verifi_time_exact = [exact_time2; exact_time5];
% falsify_time = [falsify_time2; falsify_time5];
% safe = [sum(rob_res2(:,1)==1); sum(rob_res5(:,1) == 1)];
% unsafe = [sum(rob_res2(:,1) == 0); sum(rob_res5(:,1) == 0)];
% unknown = [sum(rob_res2(:,1) == 2); sum(rob_res5(:,1) == 2)];
% 
% disp("Updated table after verification with exact analysis...");
% 
% T = table(epsilon, safe, unsafe, unknown, verify_time_approx, falsify_time, verify_time_exact);
% disp(T);
