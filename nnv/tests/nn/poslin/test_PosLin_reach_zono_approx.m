
lb = [-0.5; -0.5];
ub = [0.5; 0.5];

B = Box(lb, ub);
I1 = B.toZono;

A = [0.5 1; 1.5 -2];
I = I1.affineMap(A, []);

Z = PosLin.reach_zono_approx(I); % over-approximation using zonotope is very conservative
Z1 = PosLin.reach_zono_approx2(I); % over-approximation using new zonotope method

I2 = I.toStar; % over-approximation using star is less conservative than zonotope
S = PosLin.reach_star_approx(I2);
S2 = PosLin.reach_star_approx_fast(I2);

figure;
subplot(1, 2, 1);
I.plot;
title('Input Set');
subplot(1, 2, 2);
Z.plot;
hold on;
Z1.plot;
hold on;
S2.plot;
hold on;
S.plot;
title('Output set, inside is Star, outside is Zonotope');