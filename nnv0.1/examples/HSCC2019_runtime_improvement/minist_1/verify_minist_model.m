load mnist_model.mat;

Layers = [];

for i=1:6
    W = mnist_model.W{1,i};
    b = mnist_model.b{1, i}';
    Layers = [Layers Layer(W, b, 'ReLU')];
end

F = FFNN(Layers);

lb = mnist_model.min';
ub = mnist_model.max';


lb1 = lb;
ub1 = lb;
ub1(784) = 0.2;
ub1(783) = 0.2;


% It should be noted that the number of input state variable is very large,
% The input set should be construct using V-Representation with only small
% number of veritices. 
% Then, the affine Mapping I1 = W*I + b is very fast

% If we use the H-representation in this case, the affine Mapping is very
% slow since it need to transform the H-Representation to V-representation.
% This operation is expensive in general. 


% In this specical case study, we construct the input set using
% V-representation 

V = Reduction.getVertices(lb1, ub1); % get vertices of input set
I = Polyhedron('V', V'); % use V-representation is faster for affine mapping

I1 = Partition.partition_box(Box(lb1, ub1), 2); % lazy-approximate scheme + input partition
I11 = [];

for i=1:length(I1)
    V = Reduction.getVertices(I1(i).lb, I1(i).ub);
    I11 = [I11 Polyhedron('V', V')];
end


% number of partitions = 16

exact_runtime = [];
approx_runtime = [];
approx_partition_runtime = [];
mixing_runtime = [];

% run all schemes with number of cores from 1 to 4
Nmax = 5;
for i=1:4
    
% exact range analysis
[~, t1] = F.reach(I, 'exact', i, []); % exact scheme
exact_runtime = [exact_runtime t1];

% lazy-approximate range analysis
[~, t2] = F.reach(I, 'approx', i, []); % lazy-approximate scheme
approx_runtime = [approx_runtime t2];

% lazy-approximate + input partition method for range analysis
[~, t3] = F.reach(I11, 'approx', i, []); % lazy-approximate scheme
approx_partition_runtime = [approx_partition_runtime t3];

% mixing scheme for output range analysis
[~, t4] = F.reach(I, 'mix', i, Nmax); % choose N_max = 4
mixing_runtime = [mixing_runtime t4];

end

% runtime reduction in percentage
exact_reduction = [0];
approx_reduction = [0];
approx_partition_reduction = [0];
mixing_reduction = [0];

for i=2:4
    rd1 = 100*(exact_runtime(1) - exact_runtime(i))/ exact_runtime(1);
    rd2 = 100*(approx_runtime(1) - approx_runtime(i))/ approx_runtime(1);
    rd3 = 100*(approx_partition_runtime(1) - approx_partition_runtime(i))/ approx_partition_runtime(1);
    rd4 = 100*(mixing_runtime(1) - mixing_runtime(i))/ mixing_runtime(1);
    
    exact_reduction = [exact_reduction rd1];
    approx_reduction = [approx_reduction rd2];
    approx_partition_reduction = [approx_partition_reduction rd3];
    mixing_reduction = [mixing_reduction rd4];
end


