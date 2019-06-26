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