%% load network and input set
load mnist01.mat;

%% verify the image

numCores = 8; 
N = 100; % number of Images tested maximum is 100
%
relaxFactor = [0 0.25 0.5 0.75 1]; % relax factor
%
eps = [0.05 0.1 0.15 0.2]; % epsilon for the disturbance
%
N1 = length(relaxFactor); % timeout for eps = 0.25
N2 = length(eps);

r = zeros(N1, N2); % percentage of images that are robust
VT = zeros(N1, N2); % total verification time

% use glpk for LP optimization
dis_opt = [];
lp_solver = 'glpk';
t1 = tic;
for i=1:N1
    for j=1:N2
        [IS, Labels] = getInputSet(eps(j));
        t = tic;
        [r(i,j), ~, ~, ~, ~] = net.evaluateRBN(IS(1:N), Labels(1:N), "relax-star-area", numCores, relaxFactor(i), dis_opt, lp_solver);
        VT(i, j) = toc(t);
    end
end
total_VT = toc(t1);

save mnist01_relax_star_area_vs_ERAN.mat r VT eps;
%% Plot figures

% % robustness
fig1 = figure;
rb_0 = r(1,:);
rb_025 = r(2,:);
rb_05 = r(3,:);
rb_075 = r(4,:);
rb_1 = r(5,:);
rb_DeepZ = [1 0.91 0.47 0.05]; % result from ERAN
rb_DeepPoly = [1 0.94 0.67 0.14]; % result from ERAN
plot(eps, rb_0, '--*');
hold on;
plot(eps, rb_025, '--x');
hold on;
plot(eps, rb_05, '--o');
hold on;
plot(eps, rb_075, '--s');
hold on;
plot(eps, rb_1, '--^');
hold on;
plot(eps, rb_DeepZ, '-v');
hold on;
plot(eps, rb_DeepPoly, '-+');
legend('RF = 0', 'RF = 0.25', 'RF = 0.5', 'RF = 0.75', 'RF = 1', 'DeepZ', 'DeepPoly');
xlabel('$\epsilon$', 'interpreter', 'latex');
ylabel('Robustness');
title('Relax-star-area');
ax = gca; 
ax.FontSize = 13;

% % robustness
fig2 = figure;
vt_0 = VT(1,:);
vt_025 = VT(2,:);
vt_05 = VT(3,:);
vt_075 = VT(4,:);
vt_1 = VT(5,:);
vt_DeepZ = [243.7 259.6 277.1 299]; % result from ERAN
vt_DeepPoly = [93.4 95.4 96 94.3]; % result from ERAN
plot(eps, vt_0, '--*');
hold on;
plot(eps, vt_025, '--x');
hold on;
plot(eps, vt_05, '--o');
hold on;
plot(eps, vt_075, '--s');
hold on;
plot(eps, vt_1, '--^');
hold on;
plot(eps, vt_DeepZ, '-v');
hold on;
plot(eps, vt_DeepPoly, '-+');
legend('RF = 0', 'RF = 0.25', 'RF = 0.5', 'RF = 0.75', 'RF = 1', 'DeepZ', 'DeepPoly');
xlabel('$\epsilon$', 'interpreter', 'latex');
ylabel('Verification Time (s)');
title('Relax-star-area');
ax = gca; 
ax.FontSize = 13;
