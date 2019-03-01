center = [1; 1];
V = [1 0 1; 0 1 1];
P = ExamplePoly.randHrep('d',3);
S = Star([center V], P.A, P.b);

b = S.isEmptySet;

S1 = Star([1 1], [1; -1], [1; -2]);
b1 = S1.isEmptySet;
