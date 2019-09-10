load neural_network_information_4.mat
W1 = layer_1_weight_matrix;
b1 = layer_1_bias;
W2 = layer_2_weight_matrix';
b2 = layer_2_bias;

L1 = Layer(W1, b1, 'ReLU');
L2 = Layer(W2, b2, 'ReLU');

F = FFNN([L1, L2]);

lb = [0; 0];   % input range from Sherlock
ub = [10; 10]; % input range from Sherlock

I = Box(lb, ub);

n=5000;

x1 = (ub(1) - lb(1)).*rand(n, 1) + lb(1);
x2 = (ub(2) - lb(2)).*rand(n, 1) + lb(2);

Is = [x1'; x2'];
Y = F.sample(Is);
output = Y{1, 2}; % sampled output
range1 = [min(output) max(output)];

    
[R2, t2] = F.reach(I, 'approx', 4, []); % lazy-approximate scheme takes < 1 seconds to compute
save F_approx.mat F;
F.print('F_approx.info');
range2 = [R2.lb R2.ub];

I1 = Partition.partition_box(I, 4);
[R3, t3] = F.reach(I1, 'approx', 4, []); % lazy-approximate scheme
save F_approx_partition.mat F;
F.print('F_approx_partition.info');

R31 = Box.boxHull(R3);
range3 = [R31.lb R31.ub];


% compute conservativeness
CSV1 = 0;
[CSV2, r2] = CSV.getConservativeness(range2, range1);
[CSV3, r3] = CSV.getConservativeness(range3, range1);


