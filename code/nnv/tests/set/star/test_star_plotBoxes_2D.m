V = [1 1 0; 0 1 0; 0 0 1];
C = [1 0; -1 0; 0 1; 0 -1];
d = [1; 1; 1; 1];  % -1 <= a[1] <= 1, -1 <= a[2] <= 2

S1 = Star(V, C, d);

W = [2 1 1; 1 0 2; 0 1 0];
b = [0.5; 0.5; 0];

S2 = S1.affineMap(W, b);

figure;
S1.plot;
hold on;
S2.plot;


S = [S1 S2];

figure; 

Star.plots(S);

figure;
Star.plotBoxes_2D(S, 1, 2, 'red');

figure;
Star.plotBoxes_2D_noFill(S, 1, 2, 'red');
