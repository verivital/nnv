load 1L_abalone_nets.mat;

W1 = nnetwork.W{1, 1};
b1 = nnetwork.b{1, 1};

W2 = nnetwork.W{1, 2};
b2 = nnetwork.b{1, 2};

L1 = Layer(W1, b1, 'ReLU');
L2 = Layer(W2, b2, 'ReLU');

F = FFNN([L1, L2]);

ub = nnetwork.max;

% Input set
% x[i] = 0, i = 1:6
% lb(i) <= x[i] <= ub(i), i=7, 8

n = length(ub);

lb1 = zeros(8,1);
ub1 = zeros(8,1);
ub1(7) = ub(7);
ub1(8) = ub(8);

V = Reduction.getVertices(lb1, ub1);

% Input Set
I = Polyhedron('V', V');

% exact range analysis
[R1, t1] = F.reach(I, 'exact', 2, []); % exact scheme
R11 = Reduction.hypercubeHull(R1);
range1 = [R11.lb, R11.ub];
save F_1L_exact.mat F;


% lazy-approximate range analysis
[R2, t2] = F.reach(I, 'approx', 2, []); % lazy-approximate scheme
range2 = [R2.lb R2.ub];
save F_1L_approx.mat F;



% lazy-approximate + input partition method for range analysis
I1 = Partition.partition_box(I, 2); % lazy-approximate scheme + input partition
[R3, t3] = F.reach(I1, 'approx', 2, []); % lazy-approximate scheme
R31 = Reduction.hypercubeHull(R3);
range3 = [R31.lb R31.ub];
save F_1L_approx_partition.mat F;

% mixing scheme for output range analysis
[R4, t4] = F.reach(I, 'mix', 2, 10);
R41 = Reduction.hypercubeHull(R4);
range4 = [R41.lb, R41.ub];
save F_1L_mixing.mat F;


% compute conservativeness
CSV1 = 0;
[CSV2, r2] = CSV.getConservativeness(range2, range1);
[CSV3, r3] = CSV.getConservativeness(range3, range1);
[CSV4, r4] = CSV.getConservativeness(range4, range1);

