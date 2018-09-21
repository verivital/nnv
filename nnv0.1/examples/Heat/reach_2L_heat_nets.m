load 2L_heat_nets.mat;

W1 = nnetwork.W{1, 1};
b1 = nnetwork.b{1, 1};
W2 = nnetwork.W{1, 2};
b2 = nnetwork.b{1, 2};
W3 = nnetwork.W{1, 3};
b3 = nnetwork.b{1, 3};

L1 = Layer(W1, b1, 'ReLU');
L2 = Layer(W2, b2, 'ReLU');
L3 = Layer(W3, b3, 'Linear');

F = FFNN([L1 L2 L3]);  % heat network

I = Polyhedron('lb', nnetwork.min', 'ub', nnetwork.max');   % input set

[R, t] = F.reach(I, 'exact', 4, []); % exact reach set, *** 'approx' scheme gives an error.

save F.mat F; % save the verified network
F.print('F.info'); % print all information to a file

fig = figure;
R.plot;