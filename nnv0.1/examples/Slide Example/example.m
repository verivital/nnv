
W1 = [1 -1; 0.5 2; -1 1];
b1 = [-1; 0.5; 0];

W2 = [-2 1 1; 0.5 1 1];
b2 = [-0.5; -0.5];

L1 = Layer(W1, b1, 'ReLU'); % construct first layer
L2 = Layer(W2, b2, 'Linear');   % construct second layer

F = FFNN([L1 L2]); % construct Feedforward neural network

lb = [-2; -1]; % lower-bound vector of input set
ub = [2; 2];   % upper-bound vector of input set

I = Polyhedron('lb', lb, 'ub', ub); % construct input set

[R, t] = F.reach(I, 'exact', 4, []); % compute the exact reachable set

%[R, t] = F.reach(I, 'approx', 4, []); % compute over-approximate reachable set using lazy-approximate scheme

%[R, t] = F.reach(I, 'mix', 4, 2); % compute an over-approximate reachable set using mixing scheme


% plot reachable set
fig = figure;
subplot(1, 2, 1);
I.plot;
title('Input Set', 'FontSize', 20);
xlabel('x_1', 'FontSize', 16);
ylabel('x_2', 'FontSize', 16);

subplot(1, 2, 2)
R.plot
title('Output Set', 'FontSize', 20);
xlabel('y_1', 'FontSize', 16);
ylabel('y_2', 'FontSize', 16);

% verify safety

% unsafe region: x[1] >= 5 

U = Polyhedron('A', [-1 0], 'b', [-5]);

safe = true;
for i=1:length(R)
    R1 = R(i) & U; % check intersection of the reachable set and the unsafe set
   if ~R1.isEmptySet
       safe = false; % if the intersection is not empty -> the FFNN is unsafe
   end
end

