
W1 = [1 -1; 0.5 2; -1 1];
b1 = [-1; 0.5; 0];

W2 = [-2 1 1; 0.5 1 1];
b2 = [-0.5; -0.5];

W3 = [-1 1; -0.5 1];
b3 = [0; 1];

L1 = Layer(W1, b1, 'ReLU'); % construct first layer
L2 = Layer(W2, b2, 'ReLU');   % construct second layer
L3 = Layer(W3, b3, 'ReLU'); % construct third layer


F = FFNN([L1 L2 L3]); % construct Feedforward neural network

lb = [-2; -1]; % lower-bound vector of input set
ub = [2; 2];   % upper-bound vector of input set

I = Star(lb, ub); % construct input set

[R1, t1] = F.reach(I, 'exact', 1, []); % compute exact reach set

[R2, t2] = F.reach(I, 'approx-star', 1, []); % compute an over-approximate reachable set using star

% plot reachable set
fig = figure;
subplot(1, 3, 1);
I.plot;
title('Input Set', 'FontSize', 20);
xlabel('x_1', 'FontSize', 16);
ylabel('x_2', 'FontSize', 16);

subplot(1, 3, 2)
Star.plots(R1);
title('Output Set', 'FontSize', 20);
xlabel('y_1', 'FontSize', 16);
ylabel('y_2', 'FontSize', 16);


subplot(1, 3, 3)
R2.plot
hold on;
Star.plots(R1);
title('Output Set', 'FontSize', 20);
xlabel('y_1', 'FontSize', 16);
ylabel('y_2', 'FontSize', 16);
