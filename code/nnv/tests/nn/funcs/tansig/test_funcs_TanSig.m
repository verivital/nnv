% Create random set

I = ExamplePoly.randVrep;   
I.outerApprox;
V = [0 0; 1 0; 0 1];
I = Star(V', I.A, I.b, I.Internal.lb, I.Internal.ub); % input star

B = I.getBox; 
I_zono = B.toZono;

sample_size=25;

W = [0.5 1; -1 1];
I_zono = I_zono.affineMap(W, []);
I_star = I_zono.toStar;

%% test 1: TanSig reach approx star
S = TanSig.reach(I_star, 'approx-star');

%% test 2: TanSig reach abs dom
S = TanSig.reach(I_star, 'abs-dom');

%% test 3: TanSig reach approx zono
Z = TanSig.reach(I_zono, 'approx-zono');

%% test 4: TanSig reach star aprox
S = TanSig.reach_star_approx(I_star);

%% test 5: TanSig reach star aprox zono
Z = TanSig.reach_zono_approx(I_zono);

%% test 6: TanSig reach star aprox relax
S = TanSig.reach_star_approx(I);
relaxFactor = 0.5;
S2 = TanSig.reach_star_approx(I, 'approx-star', relaxFactor);
X = I.sample(10);
Y = TanSig.evaluate(X);
