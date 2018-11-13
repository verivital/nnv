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

desired_csv = 100; % conservativeness of 15%
k_max = 1; % maximum division level, k=1 -> divide each x[i] into 2 segments, k = 2 -> futher divide 2 segments into 4 ..
numOfCores = 4; % number of cores used in computation
n_samples = 5000; % use 5000 samples to estimate the output range

[R, t] = F.reach(I, 'exact', 4, []); % exact reach set
[R1, t1] = F.reach_approx_with_CSV_guarantee(I, desired_csv, k_max, numOfCores, n_samples); % exact scheme
[csv_vec, r, computed_range, est_range] = F.estimate_CSV(I, 5000);

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
R.plot;
hold on;
plot(y(1, :), y(2, :), 'o');