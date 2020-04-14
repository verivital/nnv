

lb = [-1; -1];
ub = [1 ; 1];

I = Star(lb, ub);

W = [2 1; 1 -1];

I = I.affineMap(W, []);

figure;
I.plot;
S = PosLin.stepReach(I, 1);
figure;
Star.plots(S);