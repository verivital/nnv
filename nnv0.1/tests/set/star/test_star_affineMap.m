center = [1; 1];
V = [1 0 1; 0 1 1];
P = ExamplePoly.randHrep('d',3);

V = [center V];
C = P.A;
d = P.b;

S = Star(V, C, d);
Z = S.getZono;
S = Star(V, C, d, Z);

W = [1 -1; 1 1];
b = [0.5; 0.5];
S1 = S.affineMap(W, b);

figure;
S1.Z.plot;
hold on;
S1.plot;