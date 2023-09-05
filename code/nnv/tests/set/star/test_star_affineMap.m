center = [1; 1];
V = [1 0 1; 0 1 1];
P = ExamplePoly.randHrep('d',3);

V = [center V];
C = P.A;
d = P.b;

S = Star(V, C, d);

W = [1 -1; 1 1];
b = [0.5; 0.5];
S1 = S.affineMap(W, b);


figure;
Star.plot(S);
hold on;
Star.plot(S1);

lb = [0; -1; 0];
ub = [1; 1; 1];

S3 = Star(lb, ub);
W = rand(3,3);
S4 = S3.affineMap(W, []);


P = Polyhedron('A', C, 'b', d);