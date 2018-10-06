load neural_network_information_0.mat
W1 = layer_1_weight_matrix;
b1 = layer_1_bias;
W2 = layer_2_weight_matrix';
b2 = layer_2_bias;

L1 = Layer(W1, b1, 'ReLU');
L2 = Layer(W2, b2, 'ReLU');

F = FFNN([L1, L2]);

lb = [-0.5; 0];
ub = [0; 0.5];

I = Polyhedron('lb', lb, 'ub', ub);

[R1, t1] = F.reach(I, 'exact', 4, []);  % exact scheme
save F_exact.mat F;
F.print('F_exact.info');
    
[R2, t2] = F.reach(I, 'approx', 4, []); % lazy-approximate scheme
save F_approx.mat F;
F.print('F_approx.info');
range2 = [R2.lb R2.ub];

I1 = Partition.partition_box(I, 4);
[R3, t3] = F.reach(I1, 'approx', 4, []); % lazy-approximate scheme
save F_approx_partition.mat F;
F.print('F_approx_partition.info');

R31 = Box.boxHull(R3);
range3 = [R31.lb R31.ub];


[R3, t3] = F.reach(I, 'mixing', 4, 100); % mixing scheme
save F_mixing.mat F;
F.print('F_mixing.info');


