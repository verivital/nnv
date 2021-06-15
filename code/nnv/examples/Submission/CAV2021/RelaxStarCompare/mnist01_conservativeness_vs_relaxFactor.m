%% load network and input set
load mnist01.mat;

%% verify the image

numCores = 8; 
N = 100; % number of Images tested maximum is 100
%
poolobj = gcp('nocreate');
delete(poolobj); % reset parpool
methods = ["relax-star-random", "relax-star-area", "relax-star-range", "relax-star-bound"];
eps = [0.1 0.15]; % epsilon for the disturbance
relaxFactor = [0 0.25 0.5 0.75 1]; % relax factor
%
N1 = length(methods);
N2 = length(eps); % timeout for eps = 0.25
N3 = length(relaxFactor);

r = zeros(N1, N2, N3); % percentage of images that are robust
VT = zeros(N1, N2, N3); % total verification time

% use glpk for LP optimization
dis_opt = [];
lp_solver = 'glpk';
t1 = tic;
for i=1:N1
    for j=1:N2
        [IS, Labels] = getInputSet(eps(j));
        for k=1:N3
            t = tic;
            [r(i,j, k), ~, ~, ~, ~] = net.evaluateRBN(IS(1:N), Labels(1:N), methods(i), numCores, relaxFactor(k), dis_opt, lp_solver);
            VT(i, j, k) = toc(t);
        end
    end
end
total_VT = toc(t1);

save mnist01_conservativeness_vs_relaxFactor.mat r VT N1 N2 N3 eps methods relaxFactor;
%% Plot figures

% robustness
fig1 = figure;
rb_rand = reshape(r(1,1,:), [1, N3]);
rb_area = reshape(r(2,1,:), [1, N3]);
rb_range = reshape(r(3,1,:), [1, N3]);
rb_bound = reshape(r(4,1,:), [1, N3]);
plot(relaxFactor, rb_rand, '-*');
hold on;
plot(relaxFactor, rb_area, '-x');
hold on;
plot(relaxFactor, rb_range, '-o');
hold on;
plot(relaxFactor, rb_bound, '-s');
legend('relax-star-random', 'relax-star-area', 'relax-star-range', 'relax-star-bound');
xlabel('Relaxation Factor (RF)');
ylabel('Robustness');
str = sprintf('\\epsilon = %.2f', eps(1));
title(str)
ax = gca; 
ax.FontSize = 13;

fig2 = figure;
rb_rand = reshape(r(1,2,:), [1, N3]);
rb_area = reshape(r(2,2,:), [1, N3]);
rb_range = reshape(r(3,2,:), [1, N3]);
rb_bound = reshape(r(4,2,:), [1, N3]);
plot(relaxFactor, rb_rand, '-*');
hold on;
plot(relaxFactor, rb_area, '-x');
hold on;
plot(relaxFactor, rb_range, '-o');
hold on;
plot(relaxFactor, rb_bound, '-s');
legend('relax-star-random', 'relax-star-area', 'relax-star-range', 'relax-star-bound');
xlabel('Relaxation Factor (RF)');
ylabel('Robustness');
str = sprintf('\\epsilon = %.2f', eps(2));
title(str)
ax = gca; 
ax.FontSize = 13;

