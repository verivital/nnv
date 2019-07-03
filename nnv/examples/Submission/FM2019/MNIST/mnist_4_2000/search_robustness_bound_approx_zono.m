load 2000_model_1layer.mat;
%load 2000model.mat;
load one_digits.mat;

N = length(W);
Layers = [];
for i=1:N
    L = LayerS(W{1,i}, b{1,i}', 'poslin');
    Layers = [Layers L];
end


F = FFNNS(Layers);

% we perform the test on 100 images of each digit
% n_inputs = size(ones, 1); % number of digit 1's images in the test: 100


input_vec = one_digits(1, :)';
n  = length(input_vec);

init_dis_bound = 0.0025; % initial bound of disturbance
tol = 0.0001; % accuracy in finding the maximum bound
max_steps = 100; % maximum number of searching step

lb_allowable = zeros(n, 1);
ub_allowable = ones(n, 1);

G1 = 1;
g1 = 0.5; 
G2 = -1;
g2 = -1.5;

U1 = HalfSpace(G1, g1);
U2 = HalfSpace(G2, g2);

un_robust_reg = [U1 U2]; % unrobust region is y < 0.5 or y > 1.5

method = 'approx-zono';
n_samples = 0; % do not search for falsified inputs
n_cores = 1;

[robustness_bound, t] = F.get_robustness_bound(input_vec, init_dis_bound, tol, max_steps, lb_allowable, ub_allowable, un_robust_reg, method, n_samples, n_cores);
