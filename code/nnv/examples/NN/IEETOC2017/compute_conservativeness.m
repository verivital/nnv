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

I = Polyhedron('lb', lb, 'ub', ub);

% select option for reachability algorithm

[R_exact, ~] = F.reach(I, 'exact', 4, []); % exact reach set


[R_approx, ~] = F.reach(I, 'approx', 1, []); % over-approximate reach set


[R_mix, ~] = F.reach(I, 'mix', 4, 800); % mixing scheme - over-approximate reach set



n = length(R_exact);
V_exact = 0;
for i=1:n
    V_exact = V_exact + R_exact(i).volume;
end
