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

N_test = 1000; % number of tests performed on each digit
n_attacked_pixels = 4; % number of local pixels are attached

% we perform the test on 100 images of each digit
n_inputs = size(ones, 1); % number of digit 1's images in the test: 100

bound = 0.05; % bound of disturbance

input_vec = ones(1, :)';
lb = input_vec;
ub = input_vec; 

for i=1:n_attacked_pixels
    attacked_pos = randi(784);% attack position is between 1-784
    if lb(attacked_pos) - bound >= 0
        lb(attacked_pos) = lb(attacked_pos) - bound;
    end
    if ub(attacked_pos) + bound <= 1
        ub(attacked_pos) = ub(attacked_pos) + bound;
    end  
end

V = Reduction.getVertices(lb, ub);
I = Polyhedron('V', V');

[R, t] = F.reach(I, 'exact', 2, []);
R1 = Reduction.hypercubeHull(R);
output_range = [R1.lb R1.ub];