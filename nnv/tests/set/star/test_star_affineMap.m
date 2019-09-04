center = [1; 1];
V = [1 0 1; 0 1 1];
P = ExamplePoly.randHrep('d',3);

V = [center V];
C = P.A;
d = P.b;


W = [1 -1; 1 1];
b = [0.5; 0.5];
S1 = S.affineMap(W, b);


figure;
S.plot;
hold on;
S1.plot;