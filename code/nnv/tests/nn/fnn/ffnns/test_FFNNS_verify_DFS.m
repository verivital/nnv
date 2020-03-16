
W1 = [1 -1; 0.5 2; -1 1];
b1 = [-1; 0.5; 0];

W2 = [-2 1 1; 0.5 1 1];
b2 = [-0.5; -0.5];

L1 = LayerS(W1, b1, 'poslin'); % construct first layer
L2 = LayerS(W2, b2, 'purelin');   % construct second layer

F = FFNNS([L1 L2]); % construct Feedforward neural network

lb = [-2; -1]; % lower-bound vector of input set
ub = [2; 2];   % upper-bound vector of input set

I = Star(lb, ub); % construct input set

%[R, t] = F.reach(I, 'exact', 4, []); % compute the exact reachable set

%[R, t] = F.reach(I, 'approx', 4, []); % compute over-approximate reachable set using lazy-approximate scheme

%[R, t] = F.reach(I, 'mix', 4, 2); % compute an over-approximate reachable set using mixing scheme


% plot reachable set
%fig = figure;
%subplot(1, 2, 1);
%I.plot;
%title('Input Set', 'FontSize', 20);
%xlabel('x_1', 'FontSize', 16);
%ylabel('x_2', 'FontSize', 16);

%subplot(1, 2, 2)
%R.plot
%title('Output Set', 'FontSize', 20);
%xlabel('y_1', 'FontSize', 16);
%ylabel('y_2', 'FontSize', 16);

% verify safety

% unsafe region: x[1] >= 5 

U = HalfSpace([-1 0], -5);

[safe, CEx] = F.verify_DFS('InputSet', I, 'UnsafeRegion', U, 'NumCores', 2)
%[safe, CEx] = F.verify_DFS('InputSet', I, 'UnsafeRegion', U)

