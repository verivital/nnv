load mnist2x500.mat;
load ones.mat;

Layers = [];
for i=1:3
    W = mnist2x250.W{1, i};
    b = mnist2x250.b{1,i}';
    L = Layer(W, b, 'ReLU');
    Layers = [Layers L];
end

F = FFNN(Layers);

attacked_pos = [10 11 12 13 14 15 16];

bound = 0.6; % bound of disturbance

input_vec = ones(1, :)';
lb = input_vec;
ub = input_vec; 

for i=1:length(attacked_pos)
    
    if lb(attacked_pos(i)) - bound >= 0
        lb(attacked_pos(i)) = lb(attacked_pos(i)) - bound;
    end
    if ub(attacked_pos(i)) + bound <= 1
        ub(attacked_pos(i)) = ub(attacked_pos(i)) + bound;
    end  
end

V = Reduction.getVertices(lb, ub);
I = Polyhedron('V', V');

n_cores = [1, 2];

runtime1 = [];
output1 = [];

% lazy-approximate scheme
for i=1:length(n_cores)
    [R, t] = F.reach(I, 'approx', n_cores(i), []);
    R1 = Reduction.hypercubeHull(R);
    range = [R1.lb R1.ub];
    runtime1 = [runtime1 t];
    output1 = vertcat(output1, range);
end

% lazy-approximate scheme + input partition
I1 = Partition.partition_box(Box(lb, ub), 1);
runtime2 = [];
output2 = [];
for i=1:length(n_cores)
    [R, t] = F.reach(I1, 'approx', n_cores(i), []);
    R1 = Reduction.hypercubeHull(R);
    range = [R1.lb R1.ub];
    runtime2 = [runtime2 t];
    output2 = vertcat(output2, range);
end


reduction1 = [];
reduction2 = [];


for i=1:2
    reduction1 = [reduction1 (runtime1(1) - runtime1(i))/runtime1(1)];
end


for i=1:2
    reduction2 = [reduction2 (runtime2(1) - runtime2(i))/runtime2(1)];
end


save result.mat;