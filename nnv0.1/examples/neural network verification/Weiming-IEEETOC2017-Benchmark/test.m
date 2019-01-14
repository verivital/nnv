
load NeuralNetwork7_3.mat;
Layers = [];
n = length(b);
for i=1:n - 1
    bi = cell2mat(b(i));
    Wi = cell2mat(W(i));
    Li = Layer(Wi, bi, 'ReLU');
    Layers = [Layers Li];
end
bn = cell2mat(b(n));
Wn = cell2mat(W(n));
Ln = Layer(Wn, bn, 'Linear');

Layers = [Layers Ln];

F = FFNN(Layers);

lb = [-1; -1; -1];
ub = [1; 1; 1];

I = Box(lb, ub);

sim_range = F.estimate_ranges(I, 5000);


