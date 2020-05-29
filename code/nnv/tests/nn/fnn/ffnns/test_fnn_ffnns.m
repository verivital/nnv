%to run this as a test, use results_fnn_ffnns=runtests('test_fnn_ffnns')
%requirements: file must start or end with test
%each test starts with two percent signs followed by the name
%shared vairables must appear before first test
%variables made by a test are not available to other tests.


%shared variables

%so there are actually variables that should be shared. Unfortuantely,
%since some of the tests choose to load things which have identical
%names, it is not practical to list them here. ah well. i think there
%is not a huge efficiency loss.

%___________________________________________________________________________________________________
%tests below originally taken from test_FFNNS_falsify.m

%% test 1: FFNNS falsify
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




%___________________________________________________________________________________________________
%tests below originally taken from test_FFNNS_isRobust.m

%% test 2: FFNNS isRobust
W = [1 1; 0 1];
b = [0; 0.5];
L = LayerS(W, b, 'poslin');
Layers = [L];


F = FFNNS(Layers);


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



%___________________________________________________________________________________________________
%tests below originally taken from test_FFNNS_isSafe.m

%% test 3: FFNNS isSafe
W = [1 1; 0 1];
b = [0; 0.5];
L = LayerS(W, b, 'poslin');
Layers = [L];


F = FFNNS(Layers);


lb = [-1; -1];
ub = [1; 1];

I = Star(lb, ub); % input set

[R, ~] = F.reach(I, 'exact-star');

G = [-1 0];
g = [-1.5];

U = HalfSpace(G, g); % unsafe region

n_samples = 100;


%[safe, t, counter_inputs] = F.isSafe(I, U, 'exact-star');
[safe, t, counter_inputs] = F.isSafe(I, U, 'approx-zono', n_samples);
%[safe, t, counter_inputs] = F.isSafe(I, U, 'approx-star', n_samples);
%[safe, t, counter_inputs] = F.isSafe(I, U, 'abs-dom', n_samples);
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


%___________________________________________________________________________________________________
%tests below originally taken from test_FFNNS_reach_star.m

%% test 4: FFNNS reach star


load NeuralNetwork7_3.mat;
Layers = [];
n = length(b);
for i=1:n - 1
    bi = cell2mat(b(i));
    Wi = cell2mat(W(i));
    Li = LayerS(Wi, bi, 'satlin');
    Layers = [Layers Li];
end
bn = cell2mat(b(n));
Wn = cell2mat(W(n));
Ln = LayerS(Wn, bn, 'purelin');

Layers = [Layers Ln];

F = FFNNS(Layers);
C = [1 0 0; -1 0 0; 0 1 0; 0 -1 0; 0 0 1; 0 0 -1];
d = [1; 1; 1; 1; 1; 1];

V = [0 0 0; 1 0 0; 0 1 0; 0 0 1];
I = Star(V', C, d); % input set as a Star set

% select option for reachability algorithm

[R, t] = F.reach(I, 4); % compute reach set using stars and 4 cores
save F.mat F; % save the verified network
F.print('F.info'); % print all information to a file


% generate some input to test the output
e = 0.25;
x = [];
y = [];
for x1=-1:e:1
    for x2=-1:e:1
        for x3=-1:e:1
            xi = [x1; x2; x3];
            yi = F.evaluate(xi);
            x = [x, xi];
            y = [y, yi];
        end
    end
end

fig = figure;
Star.plots(R);
hold on;
plot(y(1, :), y(2, :), 'o');


%___________________________________________________________________________________________________
%tests below originally taken from test_FFNNS_verify_DFS.m

%% test 5: FFNNS verify DFS


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



%___________________________________________________________________________________________________
%tests below originally taken from test_FFNNS_verify_MSG.m

%% test 6: FFNNS verify MSG


load ACASXU_run2a_1_1_batch_2000.mat;
Layers = [];
n = length(b);
for i=1:n - 1
    bi = cell2mat(b(i));
    Wi = cell2mat(W(i));
    Li = LayerS(Wi, bi, 'poslin');
    Layers = [Layers Li];
end
bn = cell2mat(b(n));
Wn = cell2mat(W(n));
Ln = LayerS(Wn, bn, 'purelin');

Layers = [Layers Ln];
F = FFNNS(Layers);

% Input Constraints
% 55947.69 <= i1(\rho) <= 60760,
% -3.14 <= i2 (\theta) <= 3.14,
%-3.14 <= i3 (\shi) <= -3.14
% 1145 <= i4 (\v_own) <= 1200, 
% 0 <= i5 (\v_in) <= 60

lb = [55947.69; -3.14; -3.14; 1145; 0];
ub = [60760; 3.14; 3.14; 1200; 60];

% normalize input
for i=1:5
    lb(i) = (lb(i) - means_for_scaling(i))/range_for_scaling(i);
    ub(i) = (ub(i) - means_for_scaling(i))/range_for_scaling(i);   
end

% unsafe region before scaling
unsafe_mat = [-1 0 0 0 0];
unsafe_vec = [-3.9911];

B = Box(lb, ub);

U = HalfSpace(unsafe_mat, unsafe_vec); %unsafe region
k = 5; % depth of search tree
sens_lb = 0.2; % 20%
[safe, VT, counterExamples] = F.verify_MSG(B, 'approx-zono', 5, 0.2, U)
