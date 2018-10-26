center = [1; 1];
V = [1 0 1; 0 1 1];
P = ExamplePoly.randHrep('d',3);

S = Star([center V], P.A, P.b);
% Halfspace: x[1] >= 1;
H = [-1 0];
g = [-1];

S1 = S.intersectHalfSpace(H, g);

figure;
S.plot;
hold on;
S1.plot;