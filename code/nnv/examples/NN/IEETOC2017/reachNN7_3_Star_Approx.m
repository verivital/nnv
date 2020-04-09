load NeuralNetwork7_3.mat;
Layers = [];
n = length(b);
for i=1:n - 1
    bi = cell2mat(b(i));
    Wi = cell2mat(W(i));
    Li = LayerS(Wi, bi, 'poslin');
    Layers = [Layers Li];
end
bn = cell2mat(b(n));
Wn = cell2mat(W(n));
Ln = LayerS(Wn, bn, 'purelin');

Layers = [Layers Ln];

F = FFNNS(Layers);
C = [1 0 0; -1 0 0; 0 1 0; 0 -1 0; 0 0 1; 0 0 -1];
ub = 0.2;
d = [ub; ub; ub; ub; ub; ub];

V = [0 0 0; 1 0 0; 0 1 0; 0 0 1];
I = Star(V', C, d); % input set as a Star set

% select option for reachability algorithm

I1 = I.getZono;

[R1, t1] = F.reach(I, 'exact-star'); % exact reachable set
[R2, t2] = F.reach(I, 'approx-star'); % over-approximate reach set using stars
[R3, t3] = F.reach(I1, 'approx-zono'); % over-approximate reachable set using zonotope

% generate some input to test the output
e = 0.05;
x = [];
y = [];
for x1=-ub:e:ub
    for x2=-ub:e:ub
        for x3=-ub:e:ub
            xi = [x1; x2; x3];
            yi = F.evaluate(xi);
            x = [x, xi];
            y = [y, yi];
        end
    end
end

fig = figure;
R3.plot;
hold on;
R2.plot;
hold on;
plot(y(1, :), y(2, :), '+');
hold on;
Star.plots(R1);

