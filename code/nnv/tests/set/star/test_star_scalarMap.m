center = [1; 1];
V = [1 0 1; 0 1 1];
P = ExamplePoly.randHrep('d',3);

S = Star([center V], P.A, P.b);
S1 = S.scalarMap(0.7);

figure;
Star.plot(S);
hold on;
Star.plot(S1);