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

% exact range analysis
[R1, t1] = F.reach(I, 'exact', 4, []); % exact scheme
R11 = Reduction.hypercubeHull(R1);
range1 = [R11.lb, R11.ub];
save F_2L_exact.mat F;


% lazy-approximate range analysis
[R2, t2] = F.reach(I, 'approx', 4, []); % lazy-approximate scheme
range2 = [R2.lb R2.ub];
save F_2L_approx.mat F;



% lazy-approximate + input partition method for range analysis
I1 = Partition.partition_box(I, 2); % lazy-approximate scheme + input partition
[R3, t3] = F.reach(I1, 'approx', 4, []); % lazy-approximate scheme
R31 = Reduction.hypercubeHull(R3);
range3 = [R31.lb R31.ub];
save F_2L_approx_partition.mat F;

% mixing scheme for output range analysis
[R4, t4] = F.reach(I, 'mix', 4, 10);
R41 = Reduction.hypercubeHull(R4);
range4 = [R41.lb, R41.ub];
save F_2L_mixing.mat F;

% compute conservativeness
CSV1 = 0;
[CSV2, r2] = CSV.getConservativeness(range2, range1);
[CSV3, r3] = CSV.getConservativeness(range3, range1);
[CSV4, r4] = CSV.getConservativeness(range4, range1);





