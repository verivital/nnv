load mnist5x50.mat;
load ones.mat;

Layers = [];
for i=1:6
    W = mnist_5x50.W{1, i};
    b = mnist_5x50.b{1,i}';
    L = Layer(W, b, 'ReLU');
    Layers = [Layers L];
end

F = FFNN(Layers);

% we perform the test on 100 images of each digit
n_inputs = size(ones, 1); % number of digit 1's images in the test: 100

attacked_pos = [1 2 3 4 5 6];

bound = 0.2; % bound of disturbance

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

n_cores = [4];
runtime = [];
output = []; 

% exact-scheme
for i=1:length(n_cores)
    [R, t] = F.reach(I, 'exact', n_cores(i), []);
    R1 = Reduction.hypercubeHull(R);
    range = [R1.lb R1.ub];
    runtime = [runtime t];
    output = vertcat(output, range);
end

