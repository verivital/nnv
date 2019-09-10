load network.mat;
Layers = [];
W = network.weights;
b = network.bias;
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

% x = [error, vref, vout, icurrent]
lb = [-3; 5; 5; 0];
ub = [3; 15; 15; 15];

% normalize input

V = Reduction.getVertices(lb, ub);

I = Polyhedron('V', V');

[R, t] = F.reach(I, 'exact', 4, []); % exact reach set
save F.mat F; % save the verified network
F.print('F.info'); % print all information to a file
