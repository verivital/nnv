load ACASXU_run2a_1_1_batch_2000.mat;
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

% Input Constraints
% 250 <= i1 <= 400, 0.2 <= i2 <= 0.4, -3.141592 <= i3 <= -3.141592 + 0.005
% 100 <= i4 <= 400, 0 <= i5 <= 400

lb = [250; 0.2; -3.141592; 100; 0.0];
ub = [400; 0.4; -3.141592 + 0.005; 400; 400];

I = Polyhedron('lb', lb, 'ub', ub);

[R, t] = F.reach(I, 'exact', 4, []); % exact reach set
save F.mat F; % save the verified network
F.print('F.info'); % print all information to a file
