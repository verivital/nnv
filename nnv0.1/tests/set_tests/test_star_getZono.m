center = [1; 1];
V = [1 0 1; 0 1 1];
P = ExamplePoly.randHrep('d',3);
S = Star([center V], P.A, P.b);
Z = S.getZono();
figure;
Z.plot;
hold on;
S.plot;
