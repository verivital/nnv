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

% select option for reachability algorithm

[R1, t1] = F.reach(I, 'exact', 4, []); % exact reach set
save F1.mat F; % save the verified network
F.print('F1.info'); % print all information to a file

[R2, t2] = F.reach(I, 'approx', 4, 300); % over-approximate reach set
save F2.mat F; % save the verified network
F.print('F2.info'); % print all information to a file

[R3, t3] = F.reach(I, 'mix', 4, 300); % mixing scheme - over-approximate reach set
save F3.mat F; % save the verified network
F.print('F3.info'); % print all information to a file


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

fig = figure; 
R2.plot;
hold on;
R1.plot;

fig = figure;
R3.plot;
hold on;
R1.plot;