% Create random set

I = ExamplePoly.randVrep;   
I.outerApprox;
V = [0 0; 1 0; 0 1];
I = Star(V', I.A, I.b, I.Internal.lb, I.Internal.ub); % input star

sample_size=25;

%% test 1: SatLins Evaluate
x = [-1.5; 0.5; 2];
y = SatLins.evaluate(x);
assert(isequal(y, [-1; 0.5; 1]))


%% test 2:  SatLins Reach exact star
S = SatLins.reach(I, 'exact-star'); % exach reach set using star

%% test 3:  SatLins Reach approx star
S = SatLins.reach(I, 'approx-star'); % over-approximate reach set using star

%% test 4:  SatLins Reach abs dom
S = SatLins.reach(I, 'abs-dom'); % over-approximate reach set using abstract-domain

%% test 5: SatLins Abstract Domain
S = SatLins.reach_abstract_domain(I); % over-approximate reach set

%% test 6:  SatLins Reach
S = SatLins.reach(I); % exach reach set

%% test 8: SatLins Reach Star Approx
S = SatLins.reach_star_approx(I); % over-approximate reach set

%% test 9: SatLins Reach Zono Approx
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
Zono.plot(I_zono);
title('Input Set');
subplot(1, 2, 2);
Zono.plot(Z);
hold on;
Star.plot(S);
hold on;
plot(Y(1, :), Y(2, :), '*b');
title('Output set, inside is Star, outside is Zonotope');

%% test 10: SatLins Step Reach 
S = SatLins.stepReach(I, 1);

%% test 11: SatLins Step Reach Abstract Domain
S = SatLins.stepReachAbstractDomain(I, 1);

%% test 12: SatLins Step Reach Abstract Domain
S = SatLins.stepReach(I,1);

%% test 13: SatLins Step Reach Star Approx
S = SatLins.stepReachStarApprox(I, 1);

%% test 14: SatLins Step Reach Zono Approx

lb = [-0.5; -0.5];
ub = [0.5; 0.5];

B = Box(lb, ub);
I1 = B.toZono;

A = [2 1; 1.5 -2];
I_zono = I1.affineMap(A, []);

figure;
Zono.plot(I_zono);
Z = SatLins.stepReachZonoApprox(I_zono, 1); % over-approximation using zonotope is very conservative
I2 = I_zono.toStar; % over-approximation using star is less conservative than zonotope
S = SatLins.stepReachStarApprox(I2,1);
figure;
Zono.plot(Z);
hold on;
Star.plot(S);
