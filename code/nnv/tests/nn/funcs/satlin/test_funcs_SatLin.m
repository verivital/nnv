% Create random set

I = ExamplePoly.randVrep;   
I.outerApprox;
V = [0 0; 1 0; 0 1];
I = Star(V', I.A, I.b, I.Internal.lb, I.Internal.ub); % input star

sample_size=25;

%% test 1: SatLin Evaluate
x = [-1.5; 0.5; 2];
y = SatLin.evaluate(x);
assert(isequal(y, [0; 0.5; 1]))

%% test 2:  SatLin Reach exact star
S = SatLin.reach(I, 'exact-star'); % exach reach set using star

%% test 3:  SatLin Reach approx star
S = SatLin.reach(I, 'approx-star'); % over-approximate reach set using star

%% test 4:  SatLin Reach abs dom
S = SatLin.reach(I, 'abs-dom'); % over-approximate reach set using abstract-domain

%% test 5: SatLin Abstract Domain
S = SatLin.reach_abstract_domain(I); % over-approximate reach set

%% test 6:  SatLin Reach
S = SatLin.reach(I); % exach reach set

%% test 8: SatLin Reach Star Approx
S = SatLin.reach_star_approx(I); % over-approximate reach set

%% test 9: SatLin Reach Zono Approx
lb = [-0.5; -0.5];
ub = [0.5; 0.5];

B = Box(lb, ub);
I1 = B.toZono;

A = [2 1; 1.5 -2];
I_zono = I1.affineMap(A, []);

Z = SatLin.reach_zono_approx(I_zono); % over-approximation using zonotope is very conservative
I2 = I_zono.toStar; % over-approximation using star is less conservative than zonotope
X = I2.sample(100);
Y = SatLin.evaluate(X);
S = SatLin.reach_star_approx(I2);
figure;
subplot(1, 2, 1);
Zono.plot(I_zono);
title('Input Set');
subplot(1, 2, 2);
Zono.plot(Z);
hold on;
Star.plot(S);
hold on;
plot(Y(1, :), Y(2, :), '*b');
title('Output set, inside is Star, outside is Zonotope');

%% test 10: SatLin Step Reach 
S = SatLin.stepReach(I, 1);

%% test 15: SatLin Step Reach Zono Approx

lb = [-0.5; -0.5];
ub = [0.5; 0.5];

B = Box(lb, ub);
I1 = B.toZono;

A = [2 1; 1.5 -2];
I_zono = I1.affineMap(A, []);

figure;
Zono.plot(I_zono);
Z = SatLin.stepReachZonoApprox(I_zono, 1); % over-approximation using zonotope is very conservative
I2 = I_zono.toStar; % over-approximation using star is less conservative than zonotope
S = SatLin.stepReach(I2, 1);
figure;
Zono.plot(Z);
hold on;
Star.plots(S);
