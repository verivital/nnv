I = ExamplePoly.randVrep;   % input set
I.outerApprox;
V = [0 0; 1 0; 0 1];
I = Star(V', I.A, I.b, I.Internal.lb, I.Internal.ub);

W = [1.5 1; 0 0.5];
b = [0.5; 0.5];
I1 = I.affineMap(W, b);


%% test 1: LayerS reach poslin
L1 = LayerS(W, b, 'poslin'); % relu

figure;
Star.plot(I);

figure;
Star.plot(I1);

S = L1.reach(I, 'exact-star'); % new Layer class
figure;
Star.plots(S);

%% test 2: LayerS reach satlin
L = LayerS(W, b, 'satlin');

tic;
S = L.reach(I, 'approx-star'); % new Layer class
toc;

figure;
Star.plots(S);

tic;
S1 = L.reach(I, 'approx-star', 'parallel');
toc;
