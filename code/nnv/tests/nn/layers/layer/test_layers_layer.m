%to run this as a test, use results_fnn_layer=runtests('test_fnn_layer')
%requirements: file must start or end with test
%each test starts with two percent signs followed by the name
%shared vairables must appear before first test
%variables made by a test are not available to other tests.


%shared variables



%___________________________________________________________________________________________________
%tests below originally taken from test_Layer.m

%% test 1: Layer

I = ExamplePoly.randVrep;   % input set
W = [1.5 1; 0 0.5];
b = [0.5; 0.5];
L = Layer(W, b, 'ReLU');

V = [0 0; 1 0; 0 1];
I1 = Star(V', I.A, I.b);

figure;
I.plot;
figure;
I1.plot;

R = L.reach_exact(I, 'single');
R1 = L.reach_exact(I1, 'single');

figure;
R.plot;
figure;
Star.plots(R1);
