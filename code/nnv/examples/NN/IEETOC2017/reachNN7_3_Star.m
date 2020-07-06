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

lb = [-1; -1; -1];
ub = -lb;
I = Star(lb, ub);


% select option for reachability algorithm

[R, t] = F.reach(I, 'exact-star', 4); % exact reach set using stars
save F.mat F; % save the verified network
F.print('F.info'); % print all information to a file


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
Star.plots(R);
hold on;
plot(y(1, :), y(2, :), 'o');

y1 = F.evaluate(lb);
y2 = F.evaluate(ub);
y3 = F.evaluate((lb+ub)/2);

hold on;
plot(y1(1,:), y1(2,:), 'x', 'Color', 'r');
hold on;
plot(y2(1,:), y2(2,:), 'x', 'Color', 'r');
hold on;
plot(y3(1,:), y3(2,:), 'x', 'Color', 'r');