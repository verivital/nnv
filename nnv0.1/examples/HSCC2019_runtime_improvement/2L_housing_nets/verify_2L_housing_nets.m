load 2L_housing_nets.mat;

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


lb = nnetwork.min';
ub = nnetwork.max';

% Input set
% x[i] = 0, i = 1:10
% lb(i) <= x[i] <= ub(i), i=11, 12

n = length(lb);

Ae = eye(n);
Ae(n-1:n, :) = [];
be = zeros(n-2, 1);

A1 = zeros(2, n);
A1(1, n) = 1;
A1(2, n) = -1;

A2 = zeros(2, n);
A2(1, n-1) = 1;
A2(2, n-1) = -1;

A = vertcat(A1, A2);
b = vertcat(ub(n), -lb(n), ub(n-1), -lb(n-1));

% Input Set
I = Polyhedron('A', A, 'b', b, 'Ae', Ae, 'be', be);

exact_runtime = [];
approx_runtime = [];
approx_partition_runtime = [];
mixing_runtime = [];

I1 = Partition.partition_box(I, 2); % lazy-approximate scheme + input partition 
% number of partitions = 16


% run all schemes with number of cores from 1 to 4
for i=1:4
    
% exact range analysis
[~, t1] = F.reach(I, 'exact', i, []); % exact scheme
exact_runtime = [exact_runtime t1];

% lazy-approximate range analysis
[~, t2] = F.reach(I, 'approx', i, []); % lazy-approximate scheme
approx_runtime = [approx_runtime t2];

% lazy-approximate + input partition method for range analysis
[~, t3] = F.reach(I1, 'approx', i, []); % lazy-approximate scheme
approx_partition_runtime = [approx_partition_runtime t3];

% mixing scheme for output range analysis
[~, t4] = F.reach(I, 'mix', i, 4); % choose N_max = 4
mixing_runtime = [mixing_runtime t4];

end











