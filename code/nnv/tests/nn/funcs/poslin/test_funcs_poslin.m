%to run this as a test, use results_fnn_poslin=runtests('test_fnn_poslin')
%requirements: file must start or end with test
%each test starts with two percent signs followed by the name
%shared vairables must appear before first test
%variables made by a test are not available to other tests.


%shared variables

I = ExamplePoly.randVrep;   
V = [0 0; 1 0; 0 1];
I = Star(V', I.A, I.b); % input star
X = I.sample(100);







%___________________________________________________________________________________________________
%tests below originally taken from test_PosLin_evaluate.m

%% test 1: PosLin evaluate


x = [-1; 0.5; 2];
y = PosLin.evaluate(x);
assert(isequal(y, [0; .5; 2]));




%___________________________________________________________________________________________________
%tests below originally taken from test_PosLin_reach.m

%% test 2: PosLin reach


lb = [-0.5; -0.5];
ub = [0.5; 0.5];

B = Box(lb, ub);
I_zono = B.toZono;

A = [0.5 1; 1.5 -2];
I_zono = I_zono.affineMap(A, []);


I1 = I_zono.toStar;
X = I1.sample(100);

figure;
I_zono.plot;
hold on;
plot(X(1, :), X(2, :), 'ob'); % sampled inputs

S1 = PosLin.reach(I1, 'exact-star'); % exact reach set using star
S2 = PosLin.reach(I1, 'approx-star'); % over-approximate reach set using star
S3 = PosLin.reach(I_zono, 'approx-zono'); % over-approximate reach set using zonotope

Y = PosLin.evaluate(X);

figure;
S3.plot;
hold on;
S2.plot;
hold on;
Star.plots(S1);
hold on;
plot(Y(1, :), Y(2, :), '*'); % sampled outputs




%___________________________________________________________________________________________________
%tests below originally taken from test_PosLin_reach_abstract_domain.m

%% test 3: PosLin reach abstract domain

figure;
I.plot;
hold on;
plot(X(1, :), X(2, :), 'ob'); % sampled inputs


S = PosLin.reach_abstract_domain(I); % over-approximate reach set
S1 = PosLin.reach_star_approx(I); % exach reach set

Y = PosLin.evaluate(X);

figure;
S.plot;
hold on;
Star.plots(S1);
hold on;
plot(Y(1, :), Y(2, :), '*'); % sampled outputs



%___________________________________________________________________________________________________
%tests below originally taken from test_PosLin_reach_star_approx.m

%% test 4: PosLin reach star approx

figure;
I.plot;
hold on;
plot(X(1, :), X(2, :), 'ob'); % sampled inputs

S = PosLin.reach_star_approx(I); % over-approximate reach set
S1 = PosLin.reach(I); % exach reach set

Y = PosLin.evaluate(X);


figure;
S.plot;
hold on;
Star.plots(S1);
hold on;
plot(Y(1, :), Y(2, :), '*'); % sampled outputs




%___________________________________________________________________________________________________
%tests below originally taken from test_PosLin_reach_star_approx_fast.m

%% test 5: PosLin reach star approx fast

figure;
I.plot;
hold on;
plot(X(1, :), X(2, :), 'ob'); % sampled inputs

S1 = PosLin.reach(I); % exach reach set
tic;
S2 = PosLin.reach_star_approx(I); 
toc;
tic;
S3 = PosLin.reach_star_approx_fast(I); % fast - over-approximate reach set
toc;
Y = PosLin.evaluate(X);

figure;
S3.plot; 
hold on;
S2.plot;
hold on;
Star.plots(S1);
hold on;
plot(Y(1, :), Y(2, :), '*'); % sampled outputs



%___________________________________________________________________________________________________
%tests below originally taken from test_PosLin_reach_star_approx_vs_zono.m

%% test 6: PosLin reach star approx vs zono

W1 = [1 -1; 0.5 2; -1 1];
b1 = [-1; 0.5; 0];

W2 = [-2 1 1; 0.5 1 1];
b2 = [-0.5; -0.5];

W3 = [2 -1; 0 1];
b3 = [1; 0];

L1 = LayerS(W1, b1, 'poslin'); % construct first layer
L2 = LayerS(W2, b2, 'poslin');   % construct second layer

lb = -rand(2, 1); % lower-bound vector of input set
ub = lb + [0.5; 0.5];   % upper-bound vector of input set

I_zono = Star(lb, ub); % construct input set
I1 = I_zono.getZono;

X1 = I_zono.affineMap(W1, b1);
BX1 = X1.getBox;
BX11 = X1.getBoxFast;
Y1 = PosLin.reach_star_approx_fast(X1);
B1 = Y1.getBoxFast;
BY1 = Y1.getBox;

XZ1 = I1.affineMap(W1, b1);
Z1 = PosLin.reach_zono_approx(XZ1);
BZ1 = Z1.getBox;

X2 = Y1.affineMap(W2, b2);
BX2 = X2.getBox;
BX22 = X2.getBoxFast;
Y2 = PosLin.reach_star_approx_fast(X2);
B2 = Y2.getBoxFast;

XZ2 = Z1.affineMap(W2, b2);
BXZ2 = XZ2.getBox;
Z2 = PosLin.reach_zono_approx(XZ2);
BZ2 = Z2.getBox;


X3 = Y2.affineMap(W3, b3);
BX3 = X3.getBox;
Y3 = PosLin.reach_star_approx_fast(X3);

XZ3 = Z2.affineMap(W3, b3);
BXZ3 = XZ3.getBox;
Z3 = PosLin.reach_zono_approx(XZ3);



figure;
Z2.plot;
hold on;
Y2.plot;

figure;
Z3.plot;
hold on;
Y3.plot;


%___________________________________________________________________________________________________
%tests below originally taken from test_PosLin_reach_zono_approx.m

%% test 7: PosLin reach zono approx

lb = [-0.5; -0.5];
ub = [0.5; 0.5];

B = Box(lb, ub);
I1 = B.toZono;

A = [0.5 1; 1.5 -2];
I_zono = I1.affineMap(A, []);

Z = PosLin.reach_zono_approx(I_zono); % over-approximation using zonotope is very conservative
Z1 = PosLin.reach_zono_approx2(I_zono); % over-approximation using new zonotope method

I2 = I_zono.toStar; % over-approximation using star is less conservative than zonotope
S = PosLin.reach_star_approx(I2);
S2 = PosLin.reach_star_approx_fast(I2);

figure;
subplot(1, 2, 1);
I_zono.plot;
title('Input Set');
subplot(1, 2, 2);
Z.plot;
hold on;
Z1.plot;
hold on;
S2.plot;
hold on;
S.plot;
title('Output set, inside is Star, outside is Zonotope');



%___________________________________________________________________________________________________
%tests below originally taken from test_PosLin_stepReach.m

%% test 8: PosLin stepReach
figure;
I.plot;
S = PosLin.stepReach2(I, 1);
figure;
Star.plots(S);




%___________________________________________________________________________________________________
%tests below originally taken from test_PosLin_stepReachAbstractDomain.m

%% test 9: PosLin stepReach Abstract Domain
figure;
I.plot;
S = PosLin.stepReachAbstractDomain(I, 1);
figure;
S.plot;



%___________________________________________________________________________________________________
%tests below originally taken from test_PosLin_stepReachStarApprox.m

%% test 10: PosLin stepReach Star Approx
figure;
I.plot;
S = PosLin.stepReachStarApprox(I, 1);
figure;
S.plot;



%___________________________________________________________________________________________________
%tests below originally taken from test_PosLin_stepReachZonoApprox.m

%% test 11: PosLin stepReach Zono Approx

lb = [-0.5; -0.5];
ub = [0.5; 0.5];

B = Box(lb, ub);
I1 = B.toZono;

A = [0.5 1; 1.5 -2];
I_zono = I1.affineMap(A, []);


figure;
I_zono.plot;
Z = PosLin.stepReachZonoApprox(I_zono, 1); % over-approximation using zonotope is very conservative
I2 = I_zono.toStar; % over-approximation using star is less conservative than zonotope
S = PosLin.stepReachStarApprox(I2,1);
figure;
Z.plot;
hold on;
S.plot;
