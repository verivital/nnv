V = [1 1 0; 0 1 0; 0 0 1];
C = [1 0; -1 0; 0 1; 0 -1];
d = [1; 1; 1; 1];  % -1 <= a[1] <= 1, -1 <= a[2] <= 2

S1 = Star(V, C, d);

W = [2 1 1; 1 0 2; 0 1 0];
b = [0.5; 0.5; 0];


V2 = [1 0; 0 1; 1 1];
C2  = [1; -1];
d2 = [0.5; 0.5];

S2 = Star(V2, C2, d2);

S3 = S1.affineMap(W, b);

S12 = S1.MinkowskiSum(S2);
S13 = S1.MinkowskiSum(S3);

figure;
S12.plot;
hold on;
S1.plot;
hold on;
S2.plot;
figure;
S12.plot;