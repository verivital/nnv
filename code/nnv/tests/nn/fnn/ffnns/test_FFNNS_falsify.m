
W = [1 1; 0 1];
b = [0; 0.5];
L = LayerS(W, b, 'poslin');
Layers = [L];


F = FFNNS(Layers);

lb = [-1; -1];
ub = [1; 1];

I = Star(lb, ub);

[R, ~] = F.reach(I);

G = [-1 0];
g = [-1.5];

U = HalfSpace(G, g);

n_samples = 1000;

counter_inputs = F.falsify(I, U, n_samples);
counter_outputs = F.sample(counter_inputs);

figure;
subplot(1, 2, 1);
I.plot;
hold on; 
plot(counter_inputs(1, :), counter_inputs(2, :), 'o');
title('Input Set and counter inputs');
subplot(1, 2, 2);
Star.plots(R);
hold on;
plot(counter_outputs(1, :), counter_outputs(2, :), 'o');
title('Output set and counter outputs');
