
lb = [-0.5; -0.5];
ub = [0.5; 0.5];

B = Box(lb, ub);
I = B.toZono;

A = [0.5 1; 1.5 -2];
I = I.affineMap(A, []);


I1 = I.toStar;
X = I1.sample(100);

figure;
I.plot;
hold on;
plot(X(1, :), X(2, :), 'ob'); % sampled inputs

S1 = PosLin.reach(I1, 'exact-star'); % exact reach set using star
S2 = PosLin.reach(I1, 'approx-star'); % over-approximate reach set using star
S3 = PosLin.reach(I, 'approx-zono'); % over-approximate reach set using zonotope

Y = PosLin.evaluate(X);

figure;
Zono.plots(S3);
hold on;
S2.plot('yellow');
hold on;
Star.plots(S1);
hold on;
plot(Y(1, :), Y(2, :), '*'); % sampled outputs
