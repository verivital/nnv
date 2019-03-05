
lb = [-0.5; -0.5];
ub = [0.5; 0.5];

B = Box(lb, ub);
I1 = B.toZono;

A = [0.5 1; 1.5 -2];
I = I1.affineMap(A, []);

Z = SatLin.reach_zono_approx(I); % over-approximation using zonotope is very conservative
I2 = I.toStar; % over-approximation using star is less conservative than zonotope
X = I2.sample(100);
Y = SatLin.evaluate(X);
S = SatLin.reach_star_approx(I2);
figure;
subplot(1, 2, 1);
I.plot;
title('Input Set');
subplot(1, 2, 2);
Z.plot;
hold on;
S.plot;
hold on;
plot(Y(1, :), Y(2, :), '*b');
title('Output set, inside is Star, outside is Zonotope');