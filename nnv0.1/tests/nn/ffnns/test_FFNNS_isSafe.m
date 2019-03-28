
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

n_samples = 500;


%[safe, t, counter_inputs] = F.isSafe(I, U, 'exact-star');
[safe, t, counter_inputs] = F.isSafe(I.getZono, U, 'approx-zono', n_samples);
counter_outputs = F.sample(counter_inputs);

figure;
subplot(1, 2, 1);
I.plot;
hold on;
%Star.plots(counter_inputs);
plot(counter_inputs(1, :), counter_inputs(2, :), 'o');
title('Input Set and counter input set');
subplot(1, 2, 2);
Star.plots(R);
hold on;
%U.plot;
plot(counter_outputs(1, :), counter_outputs(2, :), 'o');

title('Output set and unsafe region');
