c = [1; 2];
V = [0; 1];

C = [1; -1];
d = [1; 0];

% S1: x = 1: 2<= y <= 3

S1 = Star([c V], C, d);

W = [2 1; 1 -1];

S2 = S1.affineMap(W, []);

S12 = S1.convexHull(S2);

S3 = S1.convexHull_with_linearTransform(W);

figure;
%S12.plot;
%hold on;
S3.plot;
hold on;
S1.plot;
hold on;
S2.plot;
