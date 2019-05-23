% test reachability analysis of linear ODE using zonotope
% this example is from the paper: Reachability of Uncertain Linear Systems
% Using Zonotopes, Antoine Girard, HSCC 2005.

A = [-1 -4; 4 -1];
B = [1; 1];
lb = [0.9; -0.1];
ub = [1.1; 0.1];
B1 = Box(lb, ub);

I = B1.toZono(); % input set: I = [0.9, 1.1] x [-0.1, 0.1]
U = Zono(0, 0.05 * eye(1)); % ||u(t)|| <= mu = 0.05

h = 0.02; % time-step for reachability analysis
N = 100; % 100 iterations
order_max = 10; % maximum number of generators = 20

sys = LinearODE(A, B, [], []);
R = sys.reachZono(I, U, h, N, order_max);

figure;
Zono.plots(R);