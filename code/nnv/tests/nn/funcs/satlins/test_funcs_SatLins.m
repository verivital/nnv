%to run this as a test, use results_fnn_SatLins=runtests('test_fnn_SatLins')
%requirements: file must start or end with test
%each test starts with two percent signs followed by the name
%shared vairables must appear before first test
%variables made by a test are not available to other tests.


%common variables
I = ExamplePoly.randVrep;   
V = [0 0; 1 0; 0 1];
I = Star(V', I.A, I.b); % input star
X = I.sample(100)


%POTENTIAL QUESTION: do we want to draw different random variables for
%each test?

%POTENTIAL QUESTION: stepReach vs satlins, keep both? many of these
%tests seem identical except for the usage of a different method

%POTENTIAL QUSETION: need to figure out exactly what asserts to do for
%all of this


%___________________________________________________________________________________________________
%tests below originally taken from test_SatLins_Evaluate.m

%% test 1: SatLins Evaluate


x = [-1.5; 0.5; 2];
y = SatLins.evaluate(x);

assert(isequal(y, [-1; 0.5; 1]))


%___________________________________________________________________________________________________
%tests below originally taken from test_SatLins_reach.m


%% test 2:  SatLins Reach

figure;
I.plot;
hold on;
plot(X(1, :), X(2, :), 'ob'); % sampled inputs

S1 = SatLins.reach(I, 'exact-star'); % exach reach set using star
S2 = SatLins.reach(I, 'approx-star'); % over-approximate reach set using star
S3 = SatLins.reach(I, 'abs-dom'); % over-approximate reach set using abstract-domain


Y = SatLins.evaluate(X);

figure;
S3.plot;
hold on;
S2.plot;
hold on;
Star.plots(S1);
plot(Y(1, :), Y(2, :), '*'); % sampled outputs

%todo: figure our what to test here



%___________________________________________________________________________________________________
%tests below originally taken from test_SatLins_reach_abstract_domain.m


%% test 3: SatLins Abstract Domain;

figure;
I.plot;
hold on;
plot(X(1, :), X(2, :), 'ob'); % sampled inputs

S = SatLins.reach_abstract_domain(I); % over-approximate reach set
S1 = SatLins.reach(I); % exach reach set

Y = SatLins.evaluate(X);

figure;
S.plot;
hold on;
Star.plots(S1);
hold on;
plot(Y(1, :), Y(2, :), '*'); % sampled outputs



%___________________________________________________________________________________________________
%tests below originally taken from test_SatLins_reach_star_approx.m

%% test 4: SatLins Reach Star Approx

figure;
I.plot;
hold on;
plot(X(1, :), X(2, :), 'ob'); % sampled inputs

S = SatLins.reach_star_approx(I); % over-approximate reach set
S1 = SatLins.reach(I); % exach reach set

Y = SatLins.evaluate(X);

figure;
S.plot;
hold on;
Star.plots(S1);
hold on;
plot(Y(1, :), Y(2, :), '*'); % sampled outputs



%___________________________________________________________________________________________________
%tests below originally taken from test_SatLins_reach_zono_approx.m

%% test 5: SatLins Reach Zono Approx




lb = [-0.5; -0.5];
ub = [0.5; 0.5];

B = Box(lb, ub);
I1 = B.toZono;

A = [2 1; 1.5 -2];
I_zono = I1.affineMap(A, []);

Z = SatLins.reach_zono_approx(I_zono); % over-approximation using zonotope is very conservative
I2 = I_zono.toStar; % over-approximation using star is less conservative than zonotope
X = I2.sample(100);
Y = SatLins.evaluate(X);
S = SatLins.reach_star_approx(I2);
figure;
subplot(1, 2, 1);
I_zono.plot;
title('Input Set');
subplot(1, 2, 2);
Z.plot;
hold on;
S.plot;
hold on;
plot(Y(1, :), Y(2, :), '*b');
title('Output set, inside is Star, outside is Zonotope');



%___________________________________________________________________________________________________
%tests below originally taken from test_SatLins_stepReach.m

%% test 6: SatLins Step Reach 

figure;
I.plot;
S = SatLins.stepReach(I, 1);
figure;
Star.plots(S);



%___________________________________________________________________________________________________
%tests below originally taken from test_SatLins_stepReachAbstractDomain.m

%% test 7: SatLins Step Reach Abstract Domain


figure;
I.plot;
S = SatLins.stepReachAbstractDomain(I, 1);
S1 = SatLins.stepReach(I,1);
figure;
S.plot; % over-approximate reach set
hold on;
Star.plots(S1); % exact reach set




%___________________________________________________________________________________________________
%tests below originally taken from test_SatLins_stepReachStarApprox.m

%% test 8: SatLins Step Reach Star Approx

figure;
I.plot;
S = SatLins.stepReachStarApprox(I, 1);
S1 = SatLins.stepReach(I,1);
figure;
S.plot; % over-approximate reach set
hold on;
Star.plots(S1); % exact reach set



%___________________________________________________________________________________________________
%tests below originally taken from test_SatLins_stepReachZonoApprox.m

%% test 9: SatLins Step Reach Zono Approx


lb = [-0.5; -0.5];
ub = [0.5; 0.5];

B = Box(lb, ub);
I1 = B.toZono;

A = [2 1; 1.5 -2];
I_zono = I1.affineMap(A, []);


figure;
I_zono.plot;
Z = SatLins.stepReachZonoApprox(I_zono, 1); % over-approximation using zonotope is very conservative
I2 = I_zono.toStar; % over-approximation using star is less conservative than zonotope
S = SatLins.stepReachStarApprox(I2,1);
figure;
Z.plot;
hold on;
S.plot;
