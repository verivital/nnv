load NeuralNetwork5_3.mat;
Layers = [];
for i=1:5
    bi = cell2mat(b(i));
    Wi = cell2mat(W(i));
    Li = Layer(Wi, bi, 'ReLU');
    Layers = [Layers Li];
end

F = FFNN(Layers);
C = [1 0 0; -1 0 0; 0 1 0; 0 -1 0; 0 0 1; 0 0 -1];
d = [1; 1; 1; 1; 1; 1];
I = Polyhedron(C, d); % input set    -1 <= x(i) <= 1, i = 1, 2, 3

[R1, rn1, t1] = F.reach(I, 'exact'); % exact reach set
%[R2, rn2, t2] = F.reach(I, 'approx'); % over-approximate reach set

fig1 = figure;
R1.plot;
%fig2 = figure;
%R2.plot