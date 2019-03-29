
W = [1 1; 0 1];
b = [0; 0.5];
L = LayerS(W, b, 'poslin');
Layers = [L];


F = FFNNS(Layers); % network

input_vec = [1; 1]; % input vector for a single points
dis_bound = 0.5; % disturbance bound

G = [-1 0];
g = [-1.5];

U = HalfSpace(G, g); % unsafe robust region

n_samples = 100; % number of samples used to find counter examples

[~, ~, counter_inputs] = F.isRobust(input_vec, dis_bound, U, 'exact-star');
figure;
subplot(1, 2, 1);
Star.plots(counter_inputs);
title('counter input set');
subplot(1, 2, 2);
Star.plots(F.outputSet);
hold on;
U.plot;
title('Output set and unsafe region');

[~, ~, counter_inputs] = F.isRobust(input_vec, dis_bound, U, 'approx-zono', n_samples);
%[~, ~, counter_inputs] = F.isRobust(input_vec, dis_bound, U, 'approx-star', n_samples);
%[~, ~, counter_inputs] = F.isRobust(input_vec, dis_bound, U, 'abs-dom', n_samples);
counter_outputs = F.sample(counter_inputs);

figure;
subplot(1, 2, 1);
plot(counter_inputs(1, :), counter_inputs(2, :), 'o');
title('counter inputs');
subplot(1, 2, 2);
U.plot;
hold on;
plot(counter_outputs(1, :), counter_outputs(2, :), 'o');
title('Counter outputs and unsafe region');
