load 2L_abalone_nets.mat;

W1 = nnetwork.W{1, 1};
b1 = nnetwork.b{1, 1};

W2 = nnetwork.W{1, 2};
b2 = nnetwork.b{1, 2};

W3 = nnetwork.W{1, 3};
b3 = nnetwork.b{1, 3};

L1 = Layer(W1, b1, 'ReLU');
L2 = Layer(W2, b2, 'ReLU');
L3 = Layer(W3, b3, 'ReLU');


F = FFNN([L1, L2, L3]);


lb = nnetwork.min;
ub = nnetwork.max;

% Input set
% x[i] = 0, i = 1:6
% lb(i) <= x[i] <= ub(i), i=7, 8

n = length(lb);

lb1 = zeros(8,1);
ub1 = zeros(8,1);
ub1(7) = ub(7);
ub1(8) = ub(8);

V = Reduction.getVertices(lb1, ub1);
% Input Set
I = Polyhedron('V', V');

desired_csv = 15; % conservativeness of 15%
k_max = 5; % maximum division level, k=1 -> divide each x[i] into 2 segments, k = 2 -> futher divide 2 segments into 4 ..
numOfCores = 2; % number of cores used in computation
n_samples = 5000; % use 5000 samples to estimate the output range


[R1, t1] = F.reach_approx_with_CSV_guarantee(I, desired_csv, k_max, numOfCores, n_samples); % exact scheme
[csv_vec, r, computed_range, est_range] = F.estimate_CSV(I, 5000);