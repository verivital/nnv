lb = [1; 1];
ub = [2; 2];

S = Star(lb, ub);

V = S.V;
C = S.C;
d = S.d;

S2 = Star(V, C, d);