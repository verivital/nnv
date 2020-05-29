%to run this as a test, use results_fnn_LogSig=runtests('test_fnn_LogSig')
%requirements: file must start or end with test
%each test starts with two percent signs followed by the name
%shared vairables must appear before first test
%variables made by a test are not available to other tests.


%shared variables


I = ExamplePoly.randVrep;   
V = [0 0; 1 0; 0 1];
I = Star(V', I.A, I.b); % input star

%___________________________________________________________________________________________________
%tests below originally taken from test.m

%% test 1: LogSig



P1 = Polyhedron('A', R1.C, 'b', R1.d);
W = R1.V(:, 2:R1.nVar + 1);
P2 = P1.affineMap(W);



%___________________________________________________________________________________________________
%tests below originally taken from test_LogSig_reach.m

%% test 2: LogSig reach

B = I.getBox; 
I_zono = B.toZono;

W = [0.5 1; -1 1];

I_zono = I_zono.affineMap(W, []);

I_star = I_zono.toStar;
S1 = LogSig.reach(I_star, 'approx-star');
S2 = LogSig.reach(I_star, 'abs-dom');

X = I_star.sample(10);
Y = LogSig.evaluate(X);

Z = LogSig.reach(I_zono, 'approx-zono');

figure;
I_zono.plot; % input set
hold on;
plot(X(1, :), X(2, :), 'ob'); % sampled inputs

figure;
Z.plot; % zonotope reach set

figure;
Star.plots(S2, 'blue'); % abs-dom
hold on;
Star.plots(S1, 'red'); % approx-star 
hold on;
plot(Y(1, :), Y(2, :), '*'); % sampled outputs




%___________________________________________________________________________________________________
%tests below originally taken from test_LogSig_reach_star_approx.m

%% test 3: LogSig reach star approx




S = LogSig.reach_star_approx(I);
S1 = LogSig.reach_star_approx_split(I);
X = I.sample(10);
Y = LogSig.evaluate(X);



figure;
I.plot; % input set
hold on;
plot(X(1, :), X(2, :), 'ob'); % sampled inputs

figure;
Star.plots(S, 'blue'); % reach set
hold on;
Star.plots(S1, 'red');
hold on;
plot(Y(1, :), Y(2, :), '*'); % sampled outputs


%___________________________________________________________________________________________________
%tests below originally taken from test_LogSig_reach_zono_approx.m

%% test 4: LogSig reach zono approx

B = I.getBox; 
I_zono = B.toZono;

W = [0.5 1; -1 1];

I_zono = I_zono.affineMap(W, []);

I_star = I_zono.toStar;
S1 = LogSig.reach_star_approx_split(I_star);
S2 = LogSig.reach_star_approx_no_split(I_star);

X = I_star.sample(10);
Y = LogSig.evaluate(X);

Z = LogSig.reach_zono_approx(I_zono);

figure;
I.plot; % input set
hold on;
Z.plot; % zonotope reach set
hold on;
S2.plot; % star reach set = zonotope reach set
hold on;
plot(Y(1, :), Y(2, :), '*'); % sampled outputs
hold on;
plot(X(1, :), X(2, :), 'ob'); % sampled inputs




%___________________________________________________________________________________________________
%tests below originally taken from test_LogSig_stepLogSig.m

%% test 5: LogSig stepLogSig


B = I.getBox;
l = B.lb;
u = B.ub;
y_l = logsig(l);
y_u = logsig(u);
dy_l = logsig('dn', l);
dy_u = logsig('dn', u);

S = LogSig.stepLogSig_Split(I, 1, l(1), u(1), y_l(1), y_u(1), dy_l(1), dy_u(1));

figure;
Star.plots(S);

S1 = LogSig.stepLogSig_NoSplit(I, 1, l(1), u(1), y_l(1), y_u(1), dy_l(1), dy_u(1));

figure;
Star.plots(S1);

figure;
Star.plots(S1);
hold on;
Star.plots(S, 'red');
