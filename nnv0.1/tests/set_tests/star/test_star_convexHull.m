c = [1; 2];
V = [0; 1];

C = [1; -1];
d = [1; 0];

% S1: x = 1: 2<= y <= 3

S1 = Star([c V], C, d);

% S2: 0 <= x = y <= 1
c = [0; 0];
V = [1; 1];
C = [1; -1];
d = [1; 0];

S2 = Star([c V], C, d);

Z1 = S1.getZono();
Z2 = S2.getZono();
Z12 = Z1.convexHull(Z2);

S12 = S1.convexHull(S2);

figure;
S1.plot;
hold on;
S2.plot;
figure;
S12.plot;

figure;
S12.plot;
hold on;
Z12.plot;


