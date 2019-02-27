W1 = [1 -1; 0.5 2; -1 1];
b1 = [-1; 0.5; 0];

W2 = [-2 1 1; 0.5 1 1];
b2 = [-0.5; -0.5];

W3 = [-1 1; -0.5 1];
b3 = [0; 1];

L1 = Layer(W1, b1, 'ReLU'); % construct first layer
L2 = Layer(W2, b2, 'ReLU');   % construct second layer
L3 = Layer(W3, b3, 'ReLU'); % construct third layer

lb = [-2; -1]; % lower-bound vector of input set
ub = [2; 2];   % upper-bound vector of input set

I = Star(lb, ub); % construct input set


R1 = L1.reach_exact(I, 'single');

Star.plots(R1);

R11 = L1.reach_approx_star(I, 'single');

figure;
R11.plot;