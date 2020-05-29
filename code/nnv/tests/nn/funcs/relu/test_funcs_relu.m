%to run this as a test, use results_fnn_relu=runtests('test_fnn_relu')
%requirements: file must start or end with test
%each test starts with two percent signs followed by the name
%shared vairables must appear before first test
%variables made by a test are not available to other tests.


%shared variables

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



%___________________________________________________________________________________________________
%tests below originally taken from test.m

%% test 1: relu

R1 = L1.reach_exact(I, 'single');

Star.plots(R1);

R11 = L1.reach_approx_star(I, 'single');

figure;
R11.plot;


%___________________________________________________________________________________________________
%tests below originally taken from test_ReLU_reach_approx.m

%% test 2: ReLU reach approx

F = FFNN([L1 L2 L3]); % construct Feedforward neural network

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



%___________________________________________________________________________________________________
%tests below originally taken from test_ReLU_stepReach_Star.m

%% test 3: ReLU stepReach Star



V = [1 1 0; 0 1 0; 0 0 1];
C = [1 0; -1 0; 0 1; 0 -1];
d = [1; 1; 1; 1];  % -1 <= a[1] <= 1, -1 <= a[2] <= 2

I1 = Star(V, C, d);
B1 = I1.getBox;
lb = B1.lb;
ub = B1.ub;

index = 2; 
R = ReLU.stepReach_Star(I1, index, lb(index), ub(index));

figure;
I1.plot;
figure;
R(1).plot;
figure;
R(2).plot;
