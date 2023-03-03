load mnist_model.mat;
load one_digits.mat;

Layers = [];
for i=1:6
    W = mnist_model.W{1, i};
    b = mnist_model.b{1,i}';
    L = LayerS(W, b, 'poslin');
    Layers = [Layers L];
end

F = FFNNS(Layers);

% we perform the test on 100 images of each digit
% n_inputs = size(ones, 1); % number of digit 1's images in the test: 100


input_vec = one_digits(1, :)';
n  = length(input_vec);

init_dis_bound = 0.0058; % initial bound of disturbance
tol = 0.0001; % accuracy in finding the maximum bound
max_steps = 20; % maximum number of searching step

lb_allowable = zeros(n, 1);
ub_allowable = ones(n, 1);

G1 = 1;
g1 = 0.5; 
G2 = -1;
g2 = -1.5;

U1 = HalfSpace(G1, g1);
U2 = HalfSpace(G2, g2);

un_robust_reg = [U1 U2]; % unrobust region is y < 0.5 or y > 1.5

method = 'exact-star';
n_samples = 0; % do not search for falsified inputs
n_cores = 6;

%[robustness_bound, ~] = F.get_robustness_bound(input_vec, init_dis_bound, tol, max_steps, lb_allowable, ub_allowable, un_robust_reg, method, n_samples, n_cores);

[robust, t, counter_inputs] = F.isRobust(input_vec, init_dis_bound, un_robust_reg, 'exact-star', lb_allowable, ub_allowable, n_samples, n_cores);
