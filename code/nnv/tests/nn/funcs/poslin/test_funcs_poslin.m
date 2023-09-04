% Create random set

I = ExamplePoly.randVrep;
I.outerApprox;
V = [0 0; 1 0; 0 1];
I = Star(V', I.A, I.b, I.Internal.lb, I.Internal.ub); % input star

lb = [-0.5; -0.5];
ub = [0.5; 0.5];

B = Box(lb, ub);
I_zono = B.toZono;
A = [0.5 1; 1.5 -2];
I_zono = I_zono.affineMap(A, []);
I_star=I_zono.toStar;

%% test 1: PosLin evaluate
x = [-1; 0.5; 2];
y = PosLin.evaluate(x);
assert(isequal(y, [0; .5; 2]));

%% test 2: PosLin reach exact star
S = PosLin.reach(I_star, 'exact-star'); % exact reach set using star

%% test 3: PosLin reach approx star
S = PosLin.reach(I_star, 'approx-star'); % over-approximate reach set using star

%% test 4: PosLin reach approx zono
Z = PosLin.reach(I_zono, 'approx-zono'); % over-approximate reach set using zonotope

%% test 5: PosLin reach abstract domain
S = PosLin.reach_abstract_domain(I); % over-approximate reach set

%% test 7: PosLin reach star approx
S = PosLin.reach_star_approx(I); % over-approximate reach set

%% test 8: PosLin reach star approx
S = PosLin.reach(I); % exach reach set

%% test 13: PosLin reach zono approx
Z = PosLin.reach_zono_approx(I_zono); % over-approximation using zonotope is very conservative

%% test 15: PosLin reach zono approx
S = PosLin.reach_star_approx(I_star);

%% test 21: PosLin stepReach Zono Approx
S = PosLin.stepReachStarApprox(I_star,1);

