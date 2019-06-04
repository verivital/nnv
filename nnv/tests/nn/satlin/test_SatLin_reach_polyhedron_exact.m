
I = ExamplePoly.randVrep;   
V = [0 0; 1 0; 0 1];
I = Star(V', I.A, I.b); % input star
X = I.sample(100);

figure;
I.plot;
hold on;
plot(X(1, :), X(2, :), 'ob'); % sampled inputs

S1 = SatLin.reach(I, 'exact-star'); % exach reach set using star
S2 = SatLin.reach(I.toPolyhedron, 'exact-polyhedron');

Y = SatLin.evaluate(X);

figure;
Star.plots(S1);
hold on;
plot(Y(1, :), Y(2, :), '*'); % sampled outputs

figure;
S2.plot;