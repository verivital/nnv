%to run this as a test, use results_layers_layerS=runtests('test_layers_layerS')
%requirements: file must start or end with test
%each test starts with two percent signs followed by the name
%shared vairables must appear before first test
%variables made by a test are not available to other tests.


%shared variables
I = ExamplePoly.randVrep;   % input set
W = [1.5 1; 0 0.5];
b = [0.5; 0.5];
L = LayerS(W, b, 'satlin');

V = [0 0; 1 0; 0 1];
I = Star(V', I.A, I.b);

I1 = I.affineMap(W, b);





%___________________________________________________________________________________________________
%tests below originally taken from test_LayerS_reach_poslin.m

%% test 1: LayerS reach poslin



L1 = Layer(W, b, 'ReLU');

figure;
I.plot;

figure;
I1.plot;


S = L.reach(I, 'exact-star'); % new Layer class
figure;
Star.plots(S);

S1 = L1.reach_exact(I, 'single');  % old Layer class
figure;
Star.plots(S1);



%___________________________________________________________________________________________________
%tests below originally taken from test_LayerS_reach_satlin.m

%% test 2: LayerS reach satlin



figure;
I.plot;

figure;
I1.plot;

tic;
S = L.reach(I, 'approx-star'); % new Layer class
toc;

figure;
Star.plots(S);

tic;
S1 = L.reach(I, 'approx-star', 'parallel');
toc;
