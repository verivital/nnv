
lb = [-0.5; -0.5];
ub = [0.5; 0.5];

B = Box(lb, ub);
I1 = B.toZono;

A = [0.5 1; 1.5 -2];
I = I1.affineMap(A, []);


figure;
I.plot;
Z = SatLin.stepReachZonoApprox(I, 1); % over-approximation using zonotope is very conservative
I2 = I.toStar; % over-approximation using star is less conservative than zonotope
S = SatLin.stepReachStarApprox(I2,1);
figure;
Z.plot;
hold on;
S.plot;