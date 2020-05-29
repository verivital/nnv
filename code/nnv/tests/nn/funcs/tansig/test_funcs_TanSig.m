%to run this as a test, use results_fnn_TanSig=runtests('test_fnn_TanSig')
%requirements: file must start or end with test
%each test starts with two percent signs followed by the name
%shared vairables must appear before first test
%variables made by a test are not available to other tests.




I = ExamplePoly.randVrep;   
V = [0 0; 1 0; 0 1];
I = Star(V', I.A, I.b); % input star


%__________________________________________________
%below taken from test_TanSig_reach.m

%% test 1: TanSig reach

B = I.getBox; 
I_zono = B.toZono;

W = [0.5 1; -1 1];

I_zono = I_zono.affineMap(W, []);

I_star = I_zono.toStar;
S1 = TanSig.reach(I_star, 'approx-star');
S2 = TanSig.reach(I_star, 'abs-dom');

X = I_star.sample(10);
Y = TanSig.evaluate(X);

Z = TanSig.reach(I_zono, 'approx-zono');

figure;
I.plot; % input set
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





%__________________________________________________
%below taken from test_TanSig_reach_star_zono_aprox.m

%% test 2: TanSig reach star zono aprox

B = I.getBox; 
I = B.toZono;

W = [0.5 1; -1 1];

I = I.affineMap(W, []);

I_star = I.toStar;
S = TanSig.reach_star_approx(I_star);
X = I_star.sample(10);
Y = TanSig.evaluate(X);

Z = TanSig.reach_zono_approx(I);

figure;
I.plot; % input set
hold on;
Z.plot; % zonotope reach set
hold on;
S.plot; % star reach set = zonotope reach set
hold on;
plot(Y(1, :), Y(2, :), '*'); % sampled outputs
hold on;
plot(X(1, :), X(2, :), 'ob'); % sampled inputs

%__________________________________________________
%below taken from test_TanSig_reach_star_aprox.m


%% test 3: TanSig reach star aprox

S = TanSig.reach_star_approx(I);
X = I.sample(10);
Y = TanSig.evaluate(X);



figure;
I.plot; % input set
hold on;
S.plot; % reach set
hold on;
plot(Y(1, :), Y(2, :), '*'); % sampled outputs
hold on;
plot(X(1, :), X(2, :), 'ob'); % sampled inputs
