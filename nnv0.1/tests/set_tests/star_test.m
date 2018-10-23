V = [1 1 0; 0 1 0; 0 0 1];
C = [1 0; -1 0; 0 1; 0 -1];
d = [1; 1; 1; 1];  % -1 <= a[1] <= 1, -1 <= a[2] <= 2

S1 = Star(V, C, d);

W = [2 1 1; 1 0 2];
b = [0.5; 0.5];

S11 = S1.affineMap(W, b);


V2 = [1 0; 0 1; 1 1];
C2  = [1; -1];
d2 = [0.5; 0.5];

S2 = Star(V2, C2, d2);

S3 = S1.MinkowskiSum(S2);

fig = figure;
S2.plot;
hold on;
S1.plot;
hold on;
S3.plot;

B2 = S2.getBox();