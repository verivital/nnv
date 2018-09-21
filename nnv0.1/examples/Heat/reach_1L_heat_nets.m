load 1L_heat_nets.mat;

W1 = nnetwork.W{1, 1};
b1 = nnetwork.b{1, 1};
W2 = nnetwork.W{1, 2};
b2 = nnetwork.b{1, 2};


L1 = Layer(W1, b1, 'ReLU');
L2 = Layer(W2, b2, 'Linear');

F = FFNN([L1 L2]);  % heat network

I = Polyhedron('lb', nnetwork.min', 'ub', nnetwork.max');   % input set

%[R, t] = F.reach(I, 'exact', 1, []); % exact reach set
%[R, t] = F.reach(I, 'approx', 1, 300); % over-approximate reach set
[R, t] = F.reach(I, 'mix', 1, 800); % mixing scheme - over-approximate reach set
save F.mat F; % save the verified network
F.print('F.info'); % print all information to a file

fig = figure;
R.plot;
