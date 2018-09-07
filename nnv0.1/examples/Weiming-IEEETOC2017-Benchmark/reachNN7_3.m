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
C = [1 0 0; -1 0 0; 0 1 0; 0 -1 0; 0 0 1; 0 0 -1];
d = [1; 1; 1; 1; 1; 1];
I = Polyhedron(C, d);

I.plot;

% select option for reachability algorithm
%[R1, rn1, t1] = F.reach(I, 'exact', 'parallel'); % exact reach set with parallel computing
%[R1, rn1, t1] = F.reach(I, 'exact', 'single'); % exact reach set with single core

%[R1, rn1, t1] = F.reach(I, 'approx-box', 'single'); % over-approximate reach set using box 
%[R1, rn1, t1] = F.reach(I, 'approx-box', 'parallel'); % over-approximate reach set using box 
%[R1, rn1, t1] = F.reach(I, 'approx-polyhedron', 'single'); % over-approximate reach set using polyhedron
% generate some input to test the output
e = 0.25;
x = [];
y = [];
for x1=-1:e:1
    for x2=-1:e:1
        for x3=-1:e:1
            xi = [x1; x2; x3];
            yi = F.evaluate(xi);
            x = [x, xi];
            y = [y, yi];
        end
    end
end

fig = figure;
R1.plot;
hold on;
plot(y(1, :), y(2, :), 'o');